#include "converge.h"

bool check_convergence(Eigen::ArrayXd est, Eigen::ArrayXd last_est,
                       double abstol, double reltol) {
    Eigen::ArrayXd diff = (est - last_est).abs();

    return (diff.maxCoeff() <= abstol) || ((diff / est).maxCoeff() <= reltol);
}

