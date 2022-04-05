#ifndef RANDOM_H
#define RANDOM_H

#include "raceproxy_types.h"

/*
 * Sample `n` iid categorical values according to `probs`
 */
uvec rcat(int n, const vec probs);

/*
 * Sample a single categorical value according to `probs`, using uniform draw `u`
 */
int rcatp(vec probs, double u);

/*
 * Sample `n` iid vectors from a Dirichlet distribution with parameter `alpha`
 */
mat rdirichlet(int n, const vec alpha);

#endif
