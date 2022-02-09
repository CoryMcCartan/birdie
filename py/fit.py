import torch
import pyro
import numpy as np
import math
from tqdm import tqdm
from pyro.infer.importance import vectorized_importance_weights, psis_diagnostic
import py.utils as utils

def run_svi(svi, scheduler, model, guide, m_args, m_kwargs, N,
            it=5000, epoch=100, n_draws=1000, 
            it_avgs=300, lr=0.01, tol_rhat=1.2,
            silent=False):
    """
    Run SVI on `model` and `guide`, with arguments in `m_args` and `m_kwargs`
    """
    loss = np.zeros(it)
    params = None
    converge_idx = None
    opt_param = None
    elbo_step_tol = 0.005
    if not silent: pbar = tqdm(total=it)
    for j in range(it):
        loss[j] = svi.step(*m_args, **m_kwargs)
        
        if j == 0:
            params = np.empty((it, utils.get_flat_params().shape[0])) # set up param tracker
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
            if not silent: pbar.update(epoch)
            
            if converge_idx is None: # haven't convered yet
                r_hat_i = np.quantile(utils.R_hat_split(params, j, 0.5), 0.95)
                if r_hat_i <= tol_rhat or it - j <= it_avgs:
                    converge_idx = j
                    if not silent: pbar.set_description("Iterate averaging")
                else:
                    if not silent: pbar.set_description("R-hat %.2f" % r_hat_i)
                    pass
                
        if converge_idx is not None: # have converged
            opt_param = utils.accuml_params(opt_param)
            if j - converge_idx >= it_avgs:
                break
    if not silent: pbar.close()
    
    opt_param["params"] = {x: opt_param["params"][x].mean(0) for x in opt_param["params"]}
    pyro.get_param_store().set_state(opt_param)
    with pyro.plate("samples", n_draws, dim=-4):
        draws = guide(*m_args, **m_kwargs)
    draws = {x: draws[x].cpu().detach() for x in draws}
        
    lw, k = diagnose(model, guide, draws, m_args, m_kwargs, n_draws, N, nesting=3)
    
    return opt_param, draws, loss[0:j], lw, k
    
    
def diagnose(model, guide, draws, args, kwargs, n_draws=1, N=1, nesting=0):
    #log_weights, _, _ = vectorized_importance_weights(model, guide, *args, **kwargs, num_samples=n_draws, max_plate_nesting=3)
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
    
