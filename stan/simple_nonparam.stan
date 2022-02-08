////////////////////////////////////////
// NONPARAMETRIC MODEL FOR RACE IMPUTATION
// Cory McCartan, Janurary 2022
////////////////////////////////////////

data {
    int N; // individuals
    int n_x; // X
    int n_r; // races
    int n_gz; // G/Z combinations

    // observed data
    int<lower=1, upper=n_x> X[N];
    int<lower=1, upper=n_gz> GZ[N];
    simplex[n_r] pr_base[N];
    simplex[n_gz] p_gz; // p(G, Z)

    // prior
    real<lower=0> n_prior_obs;
}

parameters {
    vector<lower=0>[n_x] alpha;
    simplex[n_x] p_xrgz_raw[n_gz, n_r];
}

transformed parameters {
    matrix[n_x, n_r] p_xrgz[n_gz];
    for (i in 1:n_gz) {
        for (r in 1:n_r) {
            for (x in 1:n_x) {
                p_xrgz[i][x, r] = p_xrgz_raw[i, r][x];
            }
        }
    }
}

model {
    for (i in 1:N) {
        X[i] ~ categorical(p_xrgz[GZ[i]] * pr_base[i]);
    }

    for (i in 1:n_gz) {
        for (r in 1:n_r) {
            p_xrgz_raw[i, r] ~ dirichlet(alpha);
        }
    }
    alpha ~ gamma(1.5, 1.5 / n_prior_obs);
}

generated quantities {
    matrix[n_x, n_r] p_xr = rep_matrix(0, n_x, n_r);
    for (i in 1:n_gz) {
        p_xr += p_gz[i] * p_xrgz[i];
    }
}
