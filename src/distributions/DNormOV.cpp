//#include <config.h>
#include "DNormOV.h"

#include <TruncatedNormal.h>

#include <cmath>

#include <JRmath.h>

using std::vector;

#define MU(par) (*par[0])
#define SIGMA(par) (1/sqrt(*par[1]))
#define TAU(par) (*par[1])

DNormOV::DNormOV()
    : RScalarDist("dnormOV", 2, DIST_UNBOUNDED)
{}

bool DNormOV::checkParameterValue (vector<double const *> const &par) const
{
    return (TAU(par) > 0);
}

double
DNormOV::d(double x, vector<double const *> const &par, bool give_log) const
{
    return dnorm(x, MU(par), SIGMA(par), give_log);
}

double
DNormOV::p(double q, vector<double const *> const &par, bool lower, bool give_log)
  const
{
    return pnorm(q, MU(par), SIGMA(par), lower, give_log);
}

double 
DNormOV::q(double p, vector<double const *> const &par, bool lower, bool log_p)
  const
{
    return qnorm(p, MU(par), SIGMA(par), lower, log_p);
}

double 
DNormOV::r(vector<double const *> const &par, RNG *rng) const
{
    return rnorm(MU(par), SIGMA(par), rng);
}

double DNormOV::randomSample(vector<double const *> const &par,
			   double const *lower, double const *upper,
			   RNG *rng) const
{
    double mu = MU(par);
    double sigma = SIGMA(par);
    
    if (lower && upper) {
	double left = (*lower - mu)/sigma;
	double right = (*upper - mu)/sigma;
	return mu + sigma * inormal(left, right, rng);
    }
    else if (lower) {
	double left = (*lower - mu)/sigma;
	return mu + sigma * lnormal(left, rng);
    }
    else if (upper) {
	double right = (*upper - mu)/sigma;
	return mu + sigma * rnormal(right, rng);
    }
    else {
	return rnorm(mu, sigma, rng);
    }
}
