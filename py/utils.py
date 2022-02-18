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
    
def R_hat_split(params, i, a=0.6):
    chains = np.stack(np.split(params[split_frac(i, a)], 2, 0))
    N = chains.shape[1]
    M = chains.shape[0]
    within = chains.var(1, ddof=0).mean(0)
    between = chains.mean(1).var(0, ddof=1)
    return np.sqrt(1 + between / within)

def accuml_params(old=None):
    new = pyro.get_param_store().get_state()
    for param in new["params"]:
        new["params"][param] = new["params"][param].detach()[None, ...]
        if old is not None:
            new["params"][param] = torch.cat((old["params"][param], new["params"][param]))
    return new
    
