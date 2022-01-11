#ifndef GIBBS_H
#define GIBBS_H

#include <map>
#include <vector>
#include <iterator>
#include "raceproxy_types.h"
#include "random.h"
#include <cli/progress.h>

/*
 * Sets up the lookup table with an initial set of counts.
 */
void setup_n(lookup_GZRX &n_gzrx, const mat &alpha,
             int gz, int r, int x, int n_r);

/*
 * Compute baseline race probabilities from location, covariates, and name.
 * Returns a n_r-by-N matrix containing the probabilities of each race for each
 *   individual.
 */
mat calc_baseline_prob(int n_r, int N, const uvec &S, const uvec &GW,
                       const mat &lp_sr, const mat &lp_wgr, const vec &lp_r);

#endif
