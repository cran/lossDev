#ifndef MNORM_H_
#define MNORM_H_

#include <vector>

double MNorm_logLikelihood(double const *x, unsigned int m,
			    std::vector<double const *> const &parameters,
			    std::vector<std::vector<unsigned int> > const &dims,
			    double const *lower, double const *upper);

void MNorm_randomsample(double *x, double const *mu, double const *T,
			 bool prec, int nrow, RNG *rng);

#endif /* MNORM_H_ */
