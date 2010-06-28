//#include <config.h>
#include "DTOV.h"

#include <cmath>

#include <JRmath.h>

using std::vector;

#define MU(par) (*par[0])
#define TAU(par) (*par[1])
#define DF(par) (*par[2])

DTOV::DTOV()
    : RScalarDist("dtOV", 3, DIST_UNBOUNDED)
{}

bool DTOV::checkParameterValue (vector<double const *> const &par) const
{
    return (TAU(par) > 0 && DF(par) > 0);
}

double DTOV::d(double x, vector<double const *> const &par, bool give_log) const
{
    x = (x - MU(par)) * sqrt(TAU(par));
    if (give_log) {
	return dt(x, DF(par), 1) + log(TAU(par))/2;
    }
    else {
	return dt(x, DF(par), 0) * sqrt(TAU(par));
    }
}

double DTOV::p(double x, vector<double const *> const &par, bool lower, 
	     bool use_log) const
{
    return pt((x - MU(par)) * sqrt(TAU(par)), DF(par), lower, use_log);
}

double DTOV::q(double p, vector<double const *> const &par, bool lower, 
	     bool log_p) const
{
    return MU(par) + qt(p, DF(par), lower, log_p) / sqrt(TAU(par));
}

double DTOV::r(vector<double const *> const &par, RNG *rng) const
{
    return rt(DF(par), rng) / sqrt(TAU(par)) + MU(par);
}
