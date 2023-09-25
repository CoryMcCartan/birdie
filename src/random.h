#ifndef RANDOM_H
#define RANDOM_H

#include <Rmath.h>
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
 * Sample iid vector of length `m` from a Dirichlet distribution with parameter `alpha`
 */
VectorXd rdirichlet(const VectorXd alpha, const int m);

#endif
