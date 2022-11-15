#ifndef RANDOM_H
#define RANDOM_H

#include "birdie_types.h"

/*
 * Sample `n` iid categorical values according to `probs`
 */
ArrayXi rcat(int n, const ArrayXd probs);

/*
 * Sample a single categorical value according to `probs`, using uniform draw `u`
 */
int rcatp(const ArrayXd probs, double u);

/*
 * Sample `n` iid vectors from a Dirichlet distribution with parameter `alpha`
 */
MatrixXd rdirichlet(int n, const VectorXd alpha);

#endif
