//#include <config.h>
#include "DNormOV.h"

#include <cmath>

#include <JRmath.h>

using std::vector;

#define MU(par) (*par[0])
#define SIGMA(par) (1/sqrt(*par[1]))
#define TAU(par) (*par[1])

DNormOV::DNormOV()
    : DistScalarRmath("dnormOV", 2, DIST_UNBOUNDED, true, false)
{}

bool DNormOV::checkParameterValue (vector<double const *> const &par,
				 vector<vector<unsigned int> > const &dims) 
  const
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

