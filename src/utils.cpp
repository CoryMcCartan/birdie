#include "utils.h"

// ALL EXPORTED IN HEADER

bool check_convergence(Eigen::ArrayXd est, Eigen::ArrayXd last_est,
                       double abstol, double reltol) {
    Eigen::ArrayXd diff = (est - last_est).abs();

    return (diff.maxCoeff() <= abstol) || ((diff / est).maxCoeff() <= reltol);
}

Eigen::VectorXd to_simplex(Eigen::VectorXd y) {
    int k = y.size() + 1;
    Eigen::VectorXd z(k - 1); // break proportions
    Eigen::VectorXd x(k); // output simplex

    for (int i = 0; i < k - 1; i++) {
        z[i] = 1 / (1 + std::exp(-y[i] + std::log(k - i - 1)));
    }

    x[0] = z[0];
    double accuml = 1 - x[0];
    for (int i = 1; i < k - 1; i++) {
        x[i] = accuml * z[i];
        accuml -= x[i];
    }
    x[k - 1] = accuml;

    return x;
}

Eigen::VectorXd from_simplex(Eigen::VectorXd x) {
    int k = x.size();
    Eigen::VectorXd z(k - 1); // break proportions
    Eigen::VectorXd y(k - 1); // output simplex

    z[0] = x[0];
    double accuml = 1.0 - x[0];
    for (int i = 1; i < k - 1; i++) {
        z[i] = x[i] / accuml;
        accuml -= x[i];
    }

    for (int i = 0; i < k - 1; i++) {
        y[i] = std::log(z[i]) - std::log1p(-z[i]) + std::log(k - i - 1);
    }

    return y;
}
