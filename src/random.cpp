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

MatrixXd rdirichlet(int n, const VectorXd alpha) {
    int m = alpha.size();
    MatrixXd out(n, m);
    ArrayXd sums = ArrayXd::Zero(n);
    for (int i = 0; i < m; i++) {
        out.col(i) = as<VectorXd>(rgamma(n, alpha[i]));
        sums += out.col(i).array();
    }
    for (int i = 0; i < m; i++) {
        out.col(i) = out.col(i).array() / sums;
    }
    return out;
}
