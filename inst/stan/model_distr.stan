////////////////////////////////////////
// GIBBS SAMPLER FOR RACE IMPUTATION
// Cory McCartan, Janurary 2022
////////////////////////////////////////

data {
    int N; // individuals
    int n_x; // X
    int n_r; // races
    int n_gz; // G/Z combinations
    int n_s; // names

    // observed data
    int<lower=1, upper=n_x> X[N];
    int<lower=1, upper=n_s> S[N];
    int<lower=1, upper=n_gz> GZ[N];

    // probabilities from data and Census
    matrix[n_s, n_r] lp_sr; // p(S | R)
    matrix[n_gz, n_r] lp_gzr; // p(G,Z | R)
    vector[n_r] lp_r; // p(R)
    simplex[n_gz] p_gz; // p(G)
    simplex[n_x] p_x; // p(X)

    // prior
    real<lower=0> n_prior_obs;
}

transformed data {
    simplex[n_r] pr_base[N];
    vector[n_x] alpha = p_x * n_prior_obs;

    for (i in 1:N) {
        vector[n_r] tmp = lp_sr[S[i]]' + lp_gzr[GZ[i]]' + lp_r;
        pr_base[i] = exp(tmp - log_sum_exp(tmp));
    }
}

parameters {
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
}

generated quantities {
    matrix[n_x, n_r] p_xr = rep_matrix(0, n_x, n_r);
    vector[n_r] p_r[N];
    for (i in 1:n_gz) {
        p_xr += p_gz[i] * p_xrgz[i];
    }
    for (i in 1:N) {
        vector[n_r] tmp = pr_base[i] .* p_xrgz[GZ[i]][X[i]]';
        p_r[i] = tmp / sum(tmp);
    }
}
