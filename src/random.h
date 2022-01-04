#ifndef RANDOM_H
#define RANDOM_H

#include "raceproxy_types.h"

uvec rcat(int n, const vec probs);

int rcatp(vec probs, double u);

mat rdirichlet(int n, const vec alpha);

#endif
