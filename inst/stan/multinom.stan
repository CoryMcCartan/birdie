data {
    int<lower=0> n_y;
    int<lower=0> N;
    int<lower=0> p;
    int<lower=0> n_grp;

    matrix[N, p] X;
    matrix[N, n_y] Y;
    row_vector[N] w;
    array[N] int<lower=1, upper=n_grp> grp;

    int<lower=0, upper=1> has_int;
    real<lower=0> prior_sigma;
    real<lower=0> prior_beta;
    real<lower=0> prior_int;
}

transformed data {
    vector[n_y] ones_y = rep_vector(1, n_y);
    vector[N] int_col = rep_vector(has_int, N);
}

parameters {
    row_vector[n_y] intercept;
    matrix[p, n_y] beta;
    matrix[n_grp, n_y] u;

    vector<lower=0>[n_y] sigma_grp;
    cholesky_factor_corr[n_y] L;
}

transformed parameters {
    matrix[N, n_y] lsft;
    {
        matrix[N, n_y] linpred;
        matrix[n_y, n_y] Sigma = diag_pre_multiply(sigma_grp, L);
        linpred = int_col * intercept + X * beta + (Sigma * u[grp]')';

        // manual log softmax
        lsft = linpred - rep_matrix(log(exp(linpred) * ones_y), n_y);
    }
}

model {
    target += sum(w * (Y .* lsft));

    intercept ~ normal(0, prior_int);
    to_vector(beta) ~ normal(0, prior_beta);
    to_vector(u) ~ std_normal();

    sigma_grp ~ gamma(2.0, 2.0 ./ prior_sigma);
    L ~ lkj_corr_cholesky(2.0);
}
