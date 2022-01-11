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
    matrix<lower=0>[n_s, n_r] p_sr; // p(S | R)
    matrix<lower=0>[n_gz, n_r] p_gzr; // p(G,Z | R)
    vector<lower=0>[n_r] p_r; // p(R)
    vector<lower=0>[n_gz] p_gz; // p(G)

    // prior
    real<lower=0> n_prior_obs;
}

transformed data {
    simplex[n_r] pr_base[N];
    vector[n_x] alpha = rep_vector(n_prior_obs, n_x);

    for (i in 1:N) {
        pr_base[i] = p_sr[S[i]]' .* p_gzr[GZ[i]]' .* p_r;
        pr_base[i] /= sum(pr_base[i]);
    }
}

parameters {
    simplex[n_x] p_xrgz_raw[n_gz, n_r];
}

transformed parameters {
    // prep
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
    // log posterior
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
    for (i in 1:n_gz) {
        p_xr += p_gz[i] * p_xrgz[i];
    }
}
