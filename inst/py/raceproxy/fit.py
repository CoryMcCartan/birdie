import torch
import pyro
import numpy as np
import math
from tqdm import tqdm
import pyro.poutine as poutine
from pyro.infer.enum import get_importance_trace
#from pyro.infer.importance import vectorized_importance_weights, psis_diagnostic
from pyro.infer.importance import psis_diagnostic
import raceproxy.utils as utils

def run_svi(svi, scheduler, model, guide, m_args, m_kwargs, N,
            nesting=0, it=5000, epoch=100, n_draws=1000, 
            it_avgs=300, lr=0.01, tol_rhat=1.2,
            silent=False):
    """
    Run SVI on `model` and `guide`, with arguments in `m_args` and `m_kwargs`
    """
    loss = np.zeros(it)
    params = None
    converge_idx = None
    opt_param = None
    elbo_step_tol = 0.01
    subsamp_0 = m_kwargs["subsamp"]
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
            if lr >= 0.002 and abs(split_loss[0] - split_loss[1]) < elbo_step_tol:
                elbo_step_tol *= 0.8
                lr *= scheduler.kwargs["gamma"] 
                scheduler.step()
                if m_kwargs["subsamp"] < 5000 and N > 1.1 * m_kwargs["subsamp"]:
                    m_kwargs["subsamp"] = int(m_kwargs["subsamp"] * 1.1)
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
    
    # reduce subsampling
    m_kwargs["subsamp"] = min(N, 10000)
    
    draws, log_p, log_g, tr, _ = extract_draws(
        model, guide, *m_args, **m_kwargs, 
        num_samples=n_draws, max_plate_nesting=nesting
        )
    
    # adapted from pyro.infer.importance
    log_weights = (log_p - log_g).cpu().detach()
    log_weights -= log_weights.max()
    log_weights = torch.sort(log_weights, descending=False)[0]

    cutoff_index = - int(math.ceil(min(0.2 * n_draws, 3.0 * math.sqrt(n_draws)))) - 1
    lw_cutoff = max(math.log(1.0e-15), log_weights[cutoff_index])
    lw_tail = log_weights[log_weights > lw_cutoff]

    if len(lw_tail) < 10:
        k = float('inf')
    else:
        k, _ = pyro.ops.stats.fit_generalized_pareto(lw_tail.exp() - math.exp(lw_cutoff))
    
    # reset
    m_kwargs["subsamp"] = subsamp_0
        
    return (opt_param, draws, loss[0:j], lr, log_weights.numpy(), k, 
            log_p.cpu().detach().numpy(), log_g.cpu().detach().numpy())
            

def run_mle(svi, scheduler, model, guide, m_args, m_kwargs, N,
            nesting=0, it=5000, epoch=100, n_draws=1000, 
            it_avgs=300, lr=0.01, tol_rhat=1.2,
            silent=False):
    """
    Find MLE on `model` and `guide`, with arguments in `m_args` and `m_kwargs`
    """
    loss = np.zeros(it)
    params = None
    converge_idx = None
    opt_param = None
    elbo_step_tol = 0.01
    subsamp_0 = m_kwargs["subsamp"]
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
            if lr >= 0.002 and abs(split_loss[0] - split_loss[1]) < elbo_step_tol:
                elbo_step_tol *= 0.8
                lr *= scheduler.kwargs["gamma"] 
                scheduler.step()
                if m_kwargs["subsamp"] < 5000 and N > 1.1 * m_kwargs["subsamp"]:
                    m_kwargs["subsamp"] = int(m_kwargs["subsamp"] * 1.1)
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
    
    # reduce subsampling
    m_kwargs["subsamp"] = min(N, 10000)
    
    draws, log_p, log_g, tr, _ = extract_draws(
        model, guide, *m_args, **m_kwargs, 
        num_samples=n_draws, max_plate_nesting=nesting
        )
    
    # adapted from pyro.infer.importance
    log_weights = (log_p - log_g).cpu().detach()
    log_weights -= log_weights.max()
    log_weights = torch.sort(log_weights, descending=False)[0]

    cutoff_index = - int(math.ceil(min(0.2 * n_draws, 3.0 * math.sqrt(n_draws)))) - 1
    lw_cutoff = max(math.log(1.0e-15), log_weights[cutoff_index])
    lw_tail = log_weights[log_weights > lw_cutoff]

    if len(lw_tail) < 10:
        k = float('inf')
    else:
        k, _ = pyro.ops.stats.fit_generalized_pareto(lw_tail.exp() - math.exp(lw_cutoff))
    
    # reset
    m_kwargs["subsamp"] = subsamp_0
        
    return (opt_param, draws, loss[0:j], lr, log_weights.numpy(), k, 
            log_p.cpu().detach().numpy(), log_g.cpu().detach().numpy())
            

@torch.no_grad()
def log_prob(model, draws, args, kwargs, N=1, nesting=0, guide=False):
    with pyro.plate("draws", N, dim=-nesting-1):
        cond_model = pyro.condition(model, data=draws)
        tr = poutine.trace(cond_model, graph_type="flat").get_trace(*args, **kwargs)
        tr = poutine.util.prune_subsample_sites(tr)
        if guide:
            tr.compute_score_parts()
        else:
            tr.compute_log_prob()
        tr.pack_tensors()
    wd = tr.plate_to_symbol["draws"]
    lp = 0.0
    for site in tr.nodes.values():
        if site["type"] != "sample":
            continue
        if guide and "unconstrained" not in site["name"]:
            continue
        lpp = site["packed"]["unscaled_log_prob"] # parameter log prob
        lp += torch.einsum(lpp._pyro_dims + "->" + wd, [lpp]).cpu().detach()
    return lp


@torch.no_grad()
def extract_draws(model, guide, *args, **kwargs):
    """
    Vectorized computation of importance weights for models with static structure::
        
    :param model: probabilistic model defined as a function
    :param guide: guide used for sampling defined as a function
    :param num_samples: number of samples to draw from the guide (default 1)
    :param int max_plate_nesting: Bound on max number of nested :func:`pyro.plate` contexts.
    :returns: returns a ``(num_samples,)``-shaped tensor of importance weights
        and the model and guide traces that produced them
    """
    num_samples = kwargs.pop("num_samples", 1)
    max_plate_nesting = kwargs.pop("max_plate_nesting", None)
    normalized = kwargs.pop("normalized", False)

    if max_plate_nesting is None:
        raise ValueError("must provide max_plate_nesting")
    max_plate_nesting += 1

    def vectorize(fn):
        def _fn(*args, **kwargs):
            with pyro.plate(
                "num_particles_vectorized", num_samples, dim=-max_plate_nesting
            ):
                return fn(*args, **kwargs)

        return _fn

    model_trace, guide_trace = get_importance_trace(
        "flat", max_plate_nesting, vectorize(model), vectorize(guide), args, kwargs
    )

    guide_trace.pack_tensors()
    model_trace.pack_tensors(guide_trace.plate_to_symbol)

    wd = guide_trace.plate_to_symbol["num_particles_vectorized"]
    log_p = 0.0
    log_g = 0.0
    draws = dict()
    for site in model_trace.nodes.values():
        if site["type"] != "sample":
            continue
        log_p += torch.einsum(
            site["packed"]["unscaled_log_prob"]._pyro_dims + "->" + wd,
            [site["packed"]["unscaled_log_prob"]],
        )
        if site["value"].dim() != 1:
            draws[site["name"]] = site["value"].cpu().detach()

    for site in guide_trace.nodes.values():
        if site["type"] != "sample":
            continue
        log_g += torch.einsum(
            site["packed"]["unscaled_log_prob"]._pyro_dims + "->" + wd,
            [site["packed"]["unscaled_log_prob"]],
        )

    return draws, log_p, log_g, model_trace, guide_trace
