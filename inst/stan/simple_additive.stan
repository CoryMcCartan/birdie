////////////////////////////////////////
// ADDITIVE MODEL FOR RACE IMPUTATION
// Cory McCartan, Janurary 2022
////////////////////////////////////////

data {
    int N; // individuals
    int n_x; // X
    int n_r; // races
    int n_gz; // how many G/Z variables
    int n_gz_col; // how many total G/Z variable levels

    // observed data
    int<lower=1, upper=n_x> X[N];
    matrix[N, n_gz_col] GZ;
    int<lower=1> n_nonzero; // nonzero entries in GZ
    int<lower=1, upper=n_gz> GZ_var[n_gz_col]; // match GZ cols to variables
    simplex[n_r] pr_base[N];


    real<lower=0> prior_lp_x_scale;
    real<lower=0> prior_lp_xr_scale;
    real<lower=0> prior_beta_scale;
}

transformed data {
    vector[n_x] alpha = rep_vector(prior_lp_x_scale, n_x);
    row_vector[n_gz_col] gz_avg = rep_row_vector(1.0 / N, N) * GZ;
    /*
    // sparse matrix
    vector[n_nonzero] GZ_w = csr_extract_w(GZ);
    int GZ_v[n_nonzero] = csr_extract_v(GZ);
    int GZ_u[n_nonzero+1] = csr_extract_u(GZ);
    */
}

parameters {
    simplex[n_x] lp_x_raw;
    vector[n_x-1] lp_xr_raw[n_r];
    matrix[n_x, n_gz_col] beta_raw[n_r];
    row_vector<lower=0>[n_gz] beta_scale;
}

transformed parameters {
    matrix[n_x, n_gz_col] beta[n_r];

    for (r in 1:n_r) {
        for (x in 1:n_x) {
            beta[r, x] = beta_scale[GZ_var] .* beta_raw[r, x];
        }
    }
}

model {
    matrix[n_x, N] linpred[n_r];
    vector[n_x] lp_xr[n_r];
    for (r in 1:n_r) {
        linpred[r] = beta[r] * GZ';
        //linpred[r] = to_matrix(csr_matrix_times_vector(N, n_gz_col*n_x, GZ_w,
        //    GZ_v, GZ_u, to_vector(beta[r]')), n_x, N, 0);
        lp_xr[r] = log(lp_x_raw) + append_row(0, lp_xr_raw[r]);
    }

    for (i in 1:N) {
        matrix[n_x, n_r] p_xr_i;
        for (r in 1:n_r) {
            p_xr_i[:, r] = softmax(lp_xr[r] + linpred[r, :, i]);
        }
        X[i] ~ categorical(p_xr_i * pr_base[i]);
    }

    // priors
    lp_x_raw ~ dirichlet(alpha);
    beta_scale ~ student_t(4, 0.0, prior_beta_scale);
    for (r in 1:n_r) {
        lp_xr_raw[r] ~ student_t(4, 0.0, prior_lp_xr_scale);
        for (x in 1:n_x) {
            beta_raw[r, x] ~ std_normal();
        }
    }
}

generated quantities {
    matrix[n_x, n_r] p_xr;
    for (r in 1:n_r) {
        p_xr[:, r] = softmax(log(lp_x_raw) + append_row(0, lp_xr_raw[r])
                             + beta[r] * gz_avg');
    }
}
