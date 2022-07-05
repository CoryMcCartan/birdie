import sys
import pandas as pd
import numpy as np
sys.path.append("inst/py")
import raceproxy as rp
import pickle

X = pd.read_parquet("data-raw/python_only/X.parquet")['X']
GZ_mat = pd.read_parquet("data-raw/python_only/GZ_mat.parquet").to_numpy()
GZ_var = pd.read_parquet("data-raw/python_only/GZ_var.parquet")['GZ_var']
r_probs = pd.read_parquet("data-raw/python_only/r_probs.parquet").to_numpy()
preds = pd.read_parquet("data-raw/python_only/preds.parquet")
preds = {c: preds[c].to_numpy() for c in preds.columns}

fit = rp.pyro.fit_additive(
    X, GZ_mat, GZ_var, r_probs, preds,
    2, 1, {'x': 5, 'xr': 0.75, 'beta': 1},
    it=5000, epoch=50, subsamp=2048,
    n_draws=1000, it_avgs=400,
    n_err=100, lr=0.25, tol_rhat=1.2,
    method="svi", silent=False
    )
    
print(fit["draws"]["global"].mean(0))
    
with open("data-raw/python_only/fit.pkl", 'wb') as handle:
    pickle.dump(fit, handle)
