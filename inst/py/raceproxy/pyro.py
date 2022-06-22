import copy
import torch
import pyro
import numpy as np
from scipy.special import softmax
import pyro.distributions as dist
from pyro.infer import SVI
import pyro.infer.autoguide as autoguide
import pyro.poutine as poutine
import raceproxy.utils as utils
import raceproxy.fit as fit

from torch.profiler import profile, record_function, ProfilerActivity

# MODELS -------------------------------
        
def model_additive(X, GZ, GZ_var, pr_base, N=1, n_r=5, n_x=2, n_gz=1, n_gz_var=1,
        prior_scale={"x": 5.0, "xr": 0.75, "beta": 1.0}, subsamp=1000):
    r"""
    ADDITIVE MODEL
    
    X = (N,) tensor int array of outcomes
    GZ = (N, n_gz, 1, 1) tensor int indicator matrix
    GZ_var = (n_gz,) match GZ cols to random intercepts
    pr_base = (N, 1, n_r) tensor float array of baseline race predictions
    n_prior_obs = strength of prior on p(X)
    subsamp = subsampling size
    """
    
    p_x_dist = dist.Dirichlet(torch.ones((1, 1, 1, 1, n_x)) * prior_scale['x'])
    p_x = pyro.sample("p_x", p_x_dist)
    lp_xr_raw_dist = dist.StudentT(4.0, 0.0, torch.ones((1, 1, 1, n_r, n_x-1)) * prior_scale['xr']).to_event(1)
    
    # non-centered parametrization
    beta_scale_dist = dist.HalfCauchy(torch.ones((1, 1, n_gz_var, 1, 1)) * prior_scale['beta']).to_event(1)
    beta_raw_dist = dist.Normal(torch.zeros((1, 1, n_gz, n_r, n_x)), 1.0).to_event(1)
    with pyro.plate("n_gz_var", n_gz_var, dim=-2):
        beta_scale = pyro.sample("beta_scale", beta_scale_dist)
        
    with pyro.plate("n_r", n_r, dim=-1):
        lp_xr_raw = pyro.sample("lp_xr_raw", lp_xr_raw_dist)
        with pyro.plate("n_gz", n_gz, dim=-2):
            beta_raw = pyro.sample("beta_raw", beta_raw_dist)
            
    beta = beta_scale[..., GZ_var, :, :] * beta_raw # (1, n_gz, n_r, n_x)
    zeros_shape = list(lp_xr_raw.shape)
    zeros_shape[-1] = 1
    zeros_row = torch.zeros(zeros_shape)
    lp_xr = pyro.ops.special.safe_log(p_x) + torch.cat((zeros_row, lp_xr_raw), -1) # (1, 1, n_r, n_x)
            
    normalizer = torch.nn.Softmax(-1)
    with pyro.plate("N", N, subsample_size=subsamp) as ind:
        p_xr_i = normalizer(
            torch.einsum("ijk...,lmk...->lj...", GZ[None, ind], beta) 
            + lp_xr.squeeze(-3)) # (..., ind, n_r, n_x)
        p_x_i = (pr_base[ind] @ p_xr_i).squeeze(-2)
        
        x_dist = dist.Categorical(p_x_i[:, None, None, ...], validate_args=False)
        pyro.sample("X", x_dist, obs=X[ind])
        

def fit_additive(X, GZ, GZ_var, pr_base, preds, n_x=2, n_gz_var=1,
        prior_scale={"x": 5.0, "xr": 0.75, "beta": 1.0},
        it=200, epoch=100, subsamp=1000, n_draws=1000, 
        it_avgs=300, n_err=0, lr=0.01, tol_rhat=1.2,
        silent=False):
    # convert data from R to torch
    X = torch.tensor(X, dtype=torch.int16) - 1
    GZ = torch.tensor(GZ, dtype=torch.float)[:, :, None, None]
    GZ_var = torch.tensor(GZ_var, dtype=torch.long) - 1
    pr_base = torch.tensor(pr_base, dtype=torch.float)[:, None, :]
    n_gz = GZ.shape[1]
    N = X.shape[0]
    n_r = pr_base.shape[-1]
    
    pyro.enable_validation(False)
    optimizer = pyro.optim.AdamW({'lr': lr}, {'clip_norm': 10.0})
    make_scheduler = lambda: pyro.optim.ExponentialLR({
            "optimizer": torch.optim.Adam,
            "optim_args": {"lr": lr},
            "gamma": 0.85 # how much to multiply 'lr' every call to .step()
        })
    
    guide = autoguide.AutoNormal(model_additive, init_scale=0.01)
    elbo = pyro.infer.TraceMeanField_ELBO(ignore_jit_warnings=True)
        
    # rescale so ELBO ~= 1
    model = pyro.poutine.scale(model_additive, 1.0 / N)
    guide = pyro.poutine.scale(guide, 1.0 / N)
    
    sched = make_scheduler()
    svi = SVI(model, guide, sched, loss=elbo)
    
    # setup args once since they are reused
    m_args = (X, GZ, GZ_var, pr_base)
    m_kwargs = {"N": N, "n_r": n_r, "n_x": n_x, "n_gz": n_gz, "n_gz_var": n_gz_var,
                "prior_scale": prior_scale, "subsamp": subsamp}
    
    pyro.clear_param_store()
    opt_param, draws, loss, new_lr, lw, k, log_p, log_g = fit.run_svi(
        svi, sched, model, guide, m_args, m_kwargs, N, 3, # 3 = nesting
        it, epoch, n_draws, it_avgs, lr, tol_rhat, silent
        )
        
    # Average p_xr over GZ
    p_x = draws["p_x"].numpy()
    lp_xr_raw = draws["lp_xr_raw"].numpy()
    beta_scale = draws["beta_scale"].numpy()
    beta_raw = draws["beta_raw"].numpy()
    
    beta = beta_scale[..., GZ_var.detach().numpy(), :, :] * beta_raw 
    zero_row = np.zeros((p_x.shape[0], 1, 1, n_r, 1))
    lp_xr = np.log(p_x) + np.concatenate((zero_row, lp_xr_raw), -1) 
            
    draws_out = {}
    raw_error = {}
    raw_cov = {}
    ind = np.random.choice(N, n_err, replace=False)
    pr_base = pr_base[ind].numpy()
    normalizer = torch.nn.Softmax(-1)
    for key in preds:
        GZ_mean = np.array(preds[key], dtype=np.single)[:, None, None]
        
        p_xr = np.tensordot(GZ_mean, beta, (-3, -3)).squeeze() + lp_xr.squeeze()
        p_xr = np.exp(p_xr - p_xr.max(-1, keepdims=True))
        p_xr /= p_xr.sum(-1, keepdims=True)
        draws_out[key] = p_xr.transpose((0, 2, 1))
        
        if n_err == 0: 
            continue
        
        p_xr_i = softmax(
            np.tensordot(
                GZ.numpy()[None, ind], beta, axes=(2, 2)
            ).squeeze().transpose((1, 0, 2, 3))
            + lp_xr.squeeze(-3),
            axis=-1)  # (..., ind, n_r, n_x)
        p_xr_i = np.take_along_axis(p_xr_i, X[None, ind, None, None].numpy(), -1)  # (..., ind, n_r, 1)
        theta_tilde = (p_xr_i / (pr_base @ p_xr_i)).reshape((n_draws, 1, 1, -1))
        
        cov = ((p_xr[..., None] - p_xr[..., None].mean(0)) 
                * (theta_tilde - theta_tilde.mean(0))
            ).sum(0) / (n_draws - 1)
        
        raw_error[key] = np.sqrt(n_r * N * np.mean(np.square(cov), axis=-1)).transpose((1, 0))
        raw_cov[key] = cov
    
    return {
        "loss": loss,
        "log_weights": lw,
        "pareto_k": k,
        "log_p": log_p,
        "log_g": log_g,
        "draws": draws_out,
        "raw_error": raw_error,
        "raw_cov": raw_cov,
        "err_ind": ind,
        "beta_scale": beta_scale.squeeze()
        }

