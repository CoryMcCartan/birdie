import torch
import pyro
import numpy as np
import pyro.distributions as dist
from pyro.infer import SVI
import pyro.infer.autoguide as autoguide
import pyro.poutine as poutine
from tqdm import tqdm
import py.utils as utils

# MODELS -------------------------------

def model_nonparam(X, GZ, pr_base, N=1, n_r=5, n_x=2, n_gz=1, n_prior_obs=3.0, subsamp=1000):
    r"""
    NONPARAMETRIC MODEL
    
    X = (N,) tensor int array of outcomes
    GZ = (N,) tensor int array of covariates
    pr_base = (N, 1, n_r) tensor float array of baseline race predictions
    n_prior_obs = strength of prior on p(X)
    subsamp = subsampling size
    """
    
    alpha_shape = torch.ones((1, 1, n_x)) * 1.5
    alpha_scale = torch.ones((1, 1, n_x)) * 1.5 / n_prior_obs
    alpha_dist = dist.Gamma(alpha_shape, alpha_scale).to_event(1)
    alpha = pyro.sample("alpha", alpha_dist)
    
    with pyro.plate("n_r", n_r):
        with pyro.plate("n_gz", n_gz):
            p_xrgz_raw_dist = dist.Dirichlet(alpha)
            p_xrgz_raw = pyro.sample("p_xrgz_raw", p_xrgz_raw_dist) # (n_gz, n_r, n_x)
            
    with pyro.plate("N", N, subsample_size=subsamp) as ind:
        p_x_i = (pr_base[ind] @ p_xrgz_raw[GZ[ind]]).squeeze()
        x_dist = dist.Categorical(p_x_i, validate_args=False)
        pyro.sample("X", x_dist, obs=X[ind])
        
        
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
        
        
# FITTERS -------------------------------

def fit_nonparam(X, GZ, pr_base, p_gz, n_x=2, n_prior_obs=1.5, 
        it=200, epoch=100, subsamp=1000, lr=0.01, smooth=300, tol=0.005):
    X = torch.tensor(X, dtype=torch.int16) - 1
    GZ = torch.tensor(GZ, dtype=torch.long) - 1
    pr_base = torch.tensor(pr_base, dtype=torch.float)[:, None, :]
    n_gz = p_gz.shape[0]
    N = X.shape[0]
    n_r = pr_base.shape[-1]
    
    optimizer = pyro.optim.AdamW({'lr': lr}, {'clip_norm': 2.0})
    
    if n_gz > 20:
        guide = autoguide.AutoNormal(model_nonparam)
        elbo = pyro.infer.TraceMeanField_ELBO(ignore_jit_warnings=True)
    else:
        guide = autoguide.AutoMultivariateNormal(model_nonparam)
        elbo = pyro.infer.JitTrace_ELBO(max_plate_nesting=2, ignore_jit_warnings=True)
        
    model = pyro.poutine.scale(model_nonparam, 1.0 / X.shape[0])
    guide = pyro.poutine.scale(guide, 1.0 / X.shape[0])
    svi = SVI(model, guide, optimizer, loss=elbo)
    
    pyro.clear_param_store()
    loss = np.zeros((it,))
    slope = np.ones(it // 100)
    for j in range(it):
        loss[j] = svi.step(X, GZ, pr_base, N=N, n_r=n_r, n_x=n_x, n_gz=n_gz, 
                           n_prior_obs=n_prior_obs, subsamp=subsamp)
                           
        if j > 0 and j % 100 == 0:
            recent_loss = loss[np.arange(max(0, j-smooth), j)]
            slope[j // 100] = local_slope(recent_loss)
            #scheduler.step(recent_loss.mean())
            if (slope[j // 100] <= tol):
                break
            
            print(".", end="")
            sys.stdout.flush()
    print()
    sys.stdout.flush()
            
    with pyro.plate("samples", 2000, dim=-3):
        samples = guide(X, GZ, pr_base, N=N, n_r=n_r, n_x=n_x, n_gz=n_gz, 
                        n_prior_obs=n_prior_obs, subsamp=subsamp)
    
    p_xr = samples["p_xrgz_raw"].cpu().detach().numpy()
    p_xr = np.average(p_xr, -3, weights=p_gz)
    return {
        "loss": loss[0:j], 
        "slope": slope[0:j],
        "p_xr": p_xr.transpose((0, 2, 1))
        }

def fit_additive(X, GZ, GZ_var, pr_base, n_x=2, n_gz_var=1, 
        prior_scale={"x": 5.0, "xr": 0.75, "beta": 1.0},
        it=200, epoch=100, subsamp=1000, n_draws=1000, lr=0.01, tol_rhat=1.2):
    # convert data from R to torch
    X = torch.tensor(X, dtype=torch.int16) - 1
    GZ = torch.tensor(GZ, dtype=torch.float)[:, :, None, None]
    GZ_var = torch.tensor(GZ_var, dtype=torch.long) - 1
    pr_base = torch.tensor(pr_base, dtype=torch.float)[:, None, :]
    n_gz = GZ.shape[1]
    N = X.shape[0]
    n_r = pr_base.shape[-1]
    
    optimizer = pyro.optim.AdamW({'lr': lr}, {'clip_norm': 10.0})
    scheduler = pyro.optim.ExponentialLR({
            "optimizer": torch.optim.Adam,
            "optim_args": {"lr": lr},
            "gamma": 0.8 # how much to multiply 'lr' every call to .step()
        })
    
    guide = autoguide.AutoNormal(model_additive, init_scale=0.01)
    elbo = pyro.infer.TraceMeanField_ELBO(ignore_jit_warnings=True)
    #guide = autoguide.AutoMultivariateNormal(model_additive, init_scale=0.01)
    #elbo = pyro.infer.JitTrace_ELBO(ignore_jit_warnings=True)
        
    # rescale so ELBO ~= 1
    model = pyro.poutine.scale(model_additive, 1.0 / N)
    guide = pyro.poutine.scale(guide, 1.0 / N)
    
    svi = SVI(model, guide, scheduler, loss=elbo)
    
    # setup args once since they are reused
    m_args = (X, GZ, GZ_var, pr_base)
    m_kwargs = {"N": N, "n_r": n_r, "n_x": n_x, "n_gz": n_gz, "n_gz_var": n_gz_var,
                "prior_scale": prior_scale, "subsamp": subsamp}
    
    pyro.clear_param_store()
    loss = np.zeros(it)
    r_hat = np.ones(it // epoch)
    params = None
    converge_idx = None
    draws = None
    #pbar = tqdm(total=it)
    elbo_step_tol = 0.005
    lr = lr
    it_avgs = 300
    for j in range(it):
        loss[j] = svi.step(*m_args, **m_kwargs)
        
        if j == 0:
            params = np.empty((it, utils.get_flat_params().shape[0])) # set up param tracker
            r_hat[0] = 1e6
        elif converge_idx is None: # haven't convered yet
            params[j] = utils.get_flat_params()
            
        if j > 0 and j % epoch == 0:
            # plateau detection
            split_loss = np.stack(np.split(loss[utils.split_frac(j, 0.5)], 2, 0)).mean(1)
            if lr >= 0.001 and abs(split_loss[0] - split_loss[1]) < elbo_step_tol:
                elbo_step_tol *= 0.8
                lr *= 0.8 # this is just the tracker; change above also
                scheduler.step()
                if m_kwargs["subsamp"] < 5000 and N > 1.25 * m_kwargs["subsamp"]:
                    m_kwargs["subsamp"] = int(m_kwargs["subsamp"] * 1.25)
            #pbar.update(epoch)
            
            if converge_idx is None: # haven't convered yet
                r_hat_i = np.quantile(utils.R_hat_split(params, j, 0.5), 0.95)
                if r_hat_i <= tol_rhat or it - j <= it_avgs:
                    converge_idx = j
                    #pbar.set_description("Iterate averaging")
                #else:
                    #pbar.set_description("R-hat %.2f" % r_hat_i)
                
                r_hat[j // epoch] = r_hat_i
                
        if converge_idx is not None: # have converged
            draws = utils.accuml_params(draws)
            if j - converge_idx >= it_avgs:
                break
    #pbar.close()
    
    draws["params"] = {x: draws["params"][x].mean(0) for x in draws["params"]}
    pyro.get_param_store().set_state(draws)
    with pyro.plate("samples", n_draws, dim=-4):
        draws = guide(*m_args, **m_kwargs)
    draws = {x: draws[x].cpu().detach() for x in draws}
        
    lw, k = utils.diagnose(model, guide, draws, m_args, m_kwargs, 
            n_draws, N, nesting=3)
    
        
    # Average p_xr over GZ
    p_x = draws["p_x"].numpy()
    lp_xr_raw = draws["lp_xr_raw"].numpy()
    beta_scale = draws["beta_scale"].numpy()
    beta_raw = draws["beta_raw"].numpy()
    
    beta = beta_scale[..., GZ_var.detach().numpy(), :, :] * beta_raw 
    zero_row = np.zeros((p_x.shape[0], 1, 1, n_r, 1))
    lp_xr = np.log(p_x) + np.concatenate((zero_row, lp_xr_raw), -1) 
            
    GZ_mean = GZ.detach().numpy().mean(0)
    p_xr = np.tensordot(GZ_mean, beta, (-3, -3)).squeeze() + lp_xr.squeeze()
    p_xr = np.exp(p_xr - p_xr.max(-1, keepdims=True))
    p_xr /= p_xr.sum(-1, keepdims=True)
    
    
    return {
        "loss": loss[0:j], 
        "log_weights": lw,
        "pareto_k": k,
        "r_hat": r_hat[0:converge_idx],
        "p_xr": p_xr.transpose((0, 2, 1))
        }

