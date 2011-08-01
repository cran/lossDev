//#include <config.h>
#include "DUnifOV.h"

#include <cmath>
#include <RNG.h>

using std::vector;
using std::log;

#define LOWER(par) (*par[0])
#define UPPER(par) (*par[1])

DUnifOV::DUnifOV()
    : ScalarDist("dunifOV", 2, DIST_SPECIAL)
{}

bool  DUnifOV::checkParameterValue (vector<double const *> const &par) const
{
    return (LOWER(par) < UPPER(par));
}

double DUnifOV::logDensity(double x, PDFType type,
                           vector<double const *> const &par,
                           double const *lower, double const *upper) const
{
    return log(UPPER(par) - LOWER(par));
}

double DUnifOV::randomSample(vector<double const *> const &par, 
			   double const *lower, double const *upper,
			   RNG *rng) const
{
    return LOWER(par) + rng->uniform() * (UPPER(par) - LOWER(par));
}

double DUnifOV::typicalValue(vector<double const *> const &par,
			   double const *lower, double const *upper) const
{
    return (LOWER(par) + UPPER(par))/2;
}

double DUnifOV::l(vector<double const*> const &par) const
{
    return LOWER(par);
}

double DUnifOV::u(vector<double const*> const &par) const
{
    return UPPER(par);
}

bool DUnifOV::isSupportFixed(vector<bool> const &fixmask) const
{
    return fixmask[0] && fixmask[1]; //Lower and upper bounds fixed
}
