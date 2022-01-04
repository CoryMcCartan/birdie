#ifndef GIBBS_H
#define GIBBS_H

#include "raceproxy_types.h"
#include "random.h"

mat calc_baseline_prob(int n_r, int N, const uvec S, const uvec GW,
                       const mat lp_sr, const mat lp_wgr, const vec lp_r);

#endif
