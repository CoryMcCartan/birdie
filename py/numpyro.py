import jax.random as random
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import SVI, TraceMeanField_ELBO
import numpyro.infer.autoguide as autoguide

def model_nonparam(X, GZ, pr_base, N=1, n_r=5, n_x=2, n_gz=1, n_prior_obs=3.0, subsamp=1000):
    r"""
    X = (N,) tensor int array of outcomes
    GZ = (N,) tensor int array of covariates
    pr_base = (N,n_r) tensor float array of baseline race predictions
    n_x = number of X categories
    n_gz = number of G/Z categories
    n_prior_obs = strength of prior on p(X)
    """
    
    alpha_shape = jnp.ones(n_x) * 1.5
    alpha_scale = jnp.ones(n_x) * 1.5 / n_prior_obs
    alpha_dist = dist.Gamma(alpha_shape, alpha_scale).to_event(1)
    alpha = numpyro.sample("alpha", alpha_dist)
    
    with numpyro.plate("n_r", n_r):
        with numpyro.plate("n_gz", n_gz):
            p_xrgz_raw_dist = dist.Dirichlet(jnp.broadcast_to(alpha, (n_gz, n_r, n_x)))
            p_xrgz_raw = numpyro.sample("p_xrgz_raw", p_xrgz_raw_dist)
            
    with numpyro.plate("N", N, subsample_size=subsamp) as ind:
        p_x_i = (pr_base[ind] @ p_xrgz_raw[GZ[ind]]).squeeze()
        numpyro.sample("X", dist.Categorical(p_x_i, validate_args=False), obs=X[ind])

def fit_nonparam(X, GZ, pr_base, p_gz, n_x=2, n_prior_obs=1.5, 
                 it=200, subsamp=1000, lr=0.01, smooth=300, tol=0.005):
    X = jnp.array(X, dtype=jnp.int16) - 1
    GZ = jnp.array(GZ, dtype=jnp.int32) - 1
    pr_base = jnp.expand_dims(jnp.array(pr_base, dtype=jnp.float32), 1)
    n_gz = p_gz.shape[0]
    N = X.shape[0]
    n_r = pr_base.shape[2]
    
    optimizer = numpyro.optim.ClippedAdam(step_size=lr, clip_norm=2.0)
    
    guide = autoguide.AutoNormal(model_nonparam)
    elbo = numpyro.infer.TraceMeanField_ELBO()
    #guide = autoguide.AutoMultivariateNormal(model_nonparam)
    #elbo = pyro.infer.JitTrace_ELBO(ignore_jit_warnings=True)
        
    #model = pyro.poutine.scale(model_nonparam, 1.0 / X.shape[0])
    #guide = pyro.poutine.scale(guide, 1.0 / X.shape[0])
    svi = SVI(model_nonparam, guide, optimizer, loss=elbo)
    
    res = svi.run(random.PRNGKey(0), 2000, X, GZ, pr_base, 
            N=N, n_r=n_r, n_x=n_x, n_gz=n_gz, 
            n_prior_obs=n_prior_obs, subsamp=subsamp)
            
    return res.params
            
    #with pyro.plate("samples", 2000, dim=-3):
    #    samples = guide(X, GZ, pr_base, N=N, n_r=n_r, n_x=n_x, n_gz=n_gz, 
    #                    n_prior_obs=n_prior_obs, subsamp=subsamp)
    
    #p_xr = samples["p_xrgz_raw"].cpu().detach().numpy()
    #p_xr = np.average(p_xr, 1, weights=p_gz)
    #return {
    #    "loss": loss[0:j], 
    #    "slope": slope[0:j],
    #    "p_xr": p_xr.mean(0).T,
    #    "p_xr_low": np.quantile(p_xr, 0.05, 0).T,
    #    "p_xr_high": np.quantile(p_xr, 0.95, 0).T
    #    }
