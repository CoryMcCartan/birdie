data {
    int<lower=0> n_y;
    int<lower=0> n_r;
    int<lower=0> n_x;

    vector[n_y] X[n_x, n_r]; // sums in vector format

    real<lower=0> prior_loc_alpha;
    real<lower=0> prior_shp_alpha;
}

parameters {
    vector<lower=1>[n_y] alpha[n_r];
}

model {
    real sum_alpha;
    for (k in 1:n_r) {
        target += -n_x * sum(lgamma(alpha[k]));
        sum_alpha = sum(alpha[k]);
        target +=  n_x * (lgamma(sum_alpha) - lgamma(sum_alpha + n_y));
        for (i in 1:n_x) {
            target += sum(lgamma(X[i, k] + alpha[k]));
        }

        alpha[k] ~ gamma(prior_shp_alpha, prior_shp_alpha/prior_loc_alpha);
    }
}
