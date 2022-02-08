import torch
import pyro
import numpy as np
import math
from pyro.infer.importance import vectorized_importance_weights, psis_diagnostic

def get_flat_params():
    params = pyro.get_param_store()
    loc_params = [k for k in params.keys() if ".loc" in k]
    params = np.concatenate([params[k].cpu().detach().numpy().flatten() for k in loc_params])
    return params

def split_frac(i, a):
    return np.arange(i + 1 - 2 * (math.floor(i*a) // 2), i+1)
    
def R_hat_split(params, i, a=0.75):
    chains = np.stack(np.split(params[split_frac(i, a)], 2, 0))
    N = chains.shape[1]
    M = chains.shape[0]
    within = chains.var(1, ddof=0).mean(0)
    between = chains.mean(1).var(0, ddof=1)
    return np.sqrt(1 + between / within)

def accuml_samples(new, old=None):
    for param in new:
        new[param] = new[param].cpu().detach().numpy()[None, ...]
        if old is not None:
            new[param] = np.concatenate((old[param], new[param]))
    return new

def accuml_params(old=None):
    new = pyro.get_param_store().get_state()
    for param in new["params"]:
        new["params"][param] = new["params"][param].detach()[None, ...]
        if old is not None:
            new["params"][param] = torch.cat((old["params"][param], new["params"][param]))
    return new
    

def local_slope(y):
    x = np.arange(0, y.shape[0])
    return -np.corrcoef(x, y)[0, 1]
    
def log_prob(model, draws, args, kwargs, N=1, nesting=0, guide=False):
    if not guide:
        draws = {x: draws[x].mean(0) + 0.1*(draws[x] - draws[x].mean(0)) for x in draws}
    with pyro.plate("draws", N, dim=-nesting-1):
        cond_model = pyro.condition(model, data=draws)
        tr = pyro.poutine.trace(cond_model).get_trace(*args, **kwargs)
        tr.compute_log_prob()
        tr.pack_tensors()
    wd = tr.plate_to_symbol["draws"]
    lp = 0.0
    for site in tr.nodes.values():
        if site["type"] != "sample" or site["value"].dim() == 1:
            continue
        if guide and "unconstrained" not in site["name"]:
            continue
        lpp = site["packed"]["log_prob"] # parameter log prob
        lp += torch.einsum(lpp._pyro_dims + "->" + wd, [lpp]).cpu().detach()
    return lp
    

def diagnose(model, guide, draws, args, kwargs, n_draws=1, N=1, nesting=0):
    #log_weights, _, _ = vectorized_importance_weights(model, guide, *args, **kwargs, num_samples=draws, max_plate_nesting=3)
    log_p = log_prob(model, draws, args, kwargs, n_draws, nesting)
    log_g = log_prob(guide, draws, args, kwargs, n_draws, nesting, guide=True)
    log_weights = N * (log_p - log_g)
    
    # adapted from pyro.infer.importance
    #log_weights = N * log_weights.cpu().detach()
    log_weights -= log_weights.max()
    log_weights = torch.sort(log_weights, descending=False)[0]

    cutoff_index = - int(math.ceil(min(0.2 * n_draws, 3.0 * math.sqrt(n_draws)))) - 1
    lw_cutoff = max(math.log(1.0e-15), log_weights[cutoff_index])
    lw_tail = log_weights[log_weights > lw_cutoff]

    if len(lw_tail) < 10:
        k = float('inf')
    else:
        k, _ = pyro.ops.stats.fit_generalized_pareto(lw_tail.exp() - math.exp(lw_cutoff))
        
    return log_weights.numpy(), k
