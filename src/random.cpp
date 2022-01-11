#include "random.h"

uvec rcat(int n, const vec probs) {
    int m = probs.size();
    return as<uvec>(sample(m, n, TRUE, wrap(probs)));
}

int rcatp(vec probs, double u) {
    int m = probs.size();
    int j;
    for (j = 1; j < m; j++) probs[j] += probs[j - 1];
    u *= probs[m - 1];
    for (j = 0; j < m; j++) {
        if (probs[j] >= u) break;
    }

    return j + 1;
}

mat rdirichlet(int n, const vec alpha) {
    int m = alpha.size();
    mat out(n, m);
    for (int i = 0; i < m; i++) {
        out.col(i) = as<vec>(rgamma(n, alpha[i]));
    }
    return normalise(out, 1, 1);
}

