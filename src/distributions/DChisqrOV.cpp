//#include <config.h>
#include "DChisqrOV.h"

#include <JRmath.h>

using std::vector;

#define DF(par) (*par[0])

DChisqrOV::DChisqrOV()
    : RScalarDist("dchisqrOV", 1, DIST_POSITIVE)
{}


bool 
DChisqrOV::checkParameterValue (vector<double const *> const &par) const
{
    return (DF(par) > 0);
}

double 
DChisqrOV::d(double x, vector<double const *> const &par, bool give_log) const
{
    return dchisq(x, DF(par), give_log);
}

double 
DChisqrOV::p(double q, vector<double const *> const &par, bool lower, bool log_p)
  const
{
    return pchisq(q, DF(par), lower, log_p);
}

double
DChisqrOV::q(double p, vector<double const *> const &par, bool lower, bool log_p)
const
{
    return qchisq(p, DF(par), lower, log_p);
}

double DChisqrOV::r(vector<double const *> const &par, RNG *rng) const
{
    return rchisq(DF(par), rng);
}
