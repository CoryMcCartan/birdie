#include "random.h"

ArrayXi rcat(int n, const ArrayXd probs) {
    int m = probs.size();
    return as<ArrayXi>(sample(m, n, TRUE, wrap(probs)));
}

int rcatp(ArrayXd probs, double u) {
    int m = probs.size();
    int j;
    for (j = 1; j < m; j++) probs[j] += probs[j - 1];
    u *= probs[m - 1];
    for (j = 0; j < m; j++) {
        if (probs[j] >= u) break;
    }
    if (j == m) j = m - 1; // random floating-point overflow

    return j + 1;
}

// [[Rcpp::export]]
Eigen::VectorXi mat_rcatp(Eigen::MatrixXd probs) {
    int N = probs.rows();
    VectorXd u = as<VectorXd>(runif(N));
    VectorXi R(N);

    for (int i = 0; i < N; i++) {
        R[i] = rcatp(probs.row(i), u[i]);
    }

    return R;
}

// [[Rcpp::export]]
VectorXd rdirichlet(const VectorXd alpha, const int m) {
    VectorXd out(m);
    double sum = 0.0;
    for (int i = 0; i < m; i++) {
        out[i] = R::rgamma(alpha[i], 1);
        sum += out[i];
    }
    for (int i = 0; i < m; i++) {
        out[i] = out[i] / sum;
    }
    return out;
}
