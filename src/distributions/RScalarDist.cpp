#include "RScalarDist.h"
#include <RNG.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <stdexcept>
#include <cmath>
#include <algorithm>

using std::string;
using std::vector;
using std::length_error;
using std::logic_error;
using std::log;
using std::min;
using std::max;

double RScalarDist::calPlower(double lower, 
			      vector<double const*> const &parameters) const
{
    //P(X < lower)
    if (_discrete) {
	return p(lower - 1, parameters, true, false);
    }
    else {
	return p(lower, parameters, true, false);
    }
}

double RScalarDist::calPupper(double upper,
			      vector<double const*> const &parameters) const
{
    //P(X <= upper)
    return p(upper, parameters, true, false);
}


RScalarDist::RScalarDist(string const &name, unsigned int npar, 
			 Support support, bool discrete)
  
    : ScalarDist(name, npar, support),  _support(support), _discrete(discrete),
      _npar(npar)
{
}

double 
RScalarDist::typicalValue(vector<double const *> const &parameters,
			  double const *lower, double const *upper) const
{
    double llimit = l(parameters), ulimit = u(parameters);
    double plower = 0, pupper = 1;
    
    if (lower) {
	llimit = max(llimit, *lower);
	plower = calPlower(llimit, parameters);
    }

    if (upper) {
	ulimit = min(ulimit, *upper);
	pupper = calPupper(ulimit, parameters);
    }
    
    double pmed = (plower + pupper)/2;
    double med = q(pmed, parameters, true, false);	

    //Calculate the log densities
    double dllimit = d(llimit, parameters, true);
    double dulimit = d(ulimit, parameters, true);
    double dmed = d(med, parameters, true);

    //Pick the median if it has the highest density, otherwise pick
    //a point near to (but not on) the boundary
    if (dmed >= dllimit && dmed >= dulimit) {
	return med;
    }
    else if (dulimit > dllimit) {
	return q(0.1 * plower + 0.9 * pupper, parameters, true, false);
    }
    else {
	return q(0.9 * plower + 0.1 * pupper, parameters, true, false);
    }
}

double 
RScalarDist::logLikelihood(double x, vector<double const *> const &parameters,
			   double const *lower, double const *upper) const
{
    double loglik =  d(x, parameters, true);

    if (lower || upper) {

	if (lower && x < *lower)
	    return JAGS_NEGINF;
	if (upper && x > *upper)
	    return JAGS_NEGINF;
	if (upper && lower && *upper < *lower)
	    return JAGS_NEGINF;

	//Make adjustment for discrete-valued distributions
	double ll = 0;
	if (lower) {
	    ll = _discrete ? *lower - 1 : *lower;
	}

	/* In theory, we just have to subtract log[P(lower <= X <=
           upper)] from the log likelihood. But we need to work around
           numerical problems. */

	bool have_lower = lower && p(ll, parameters, true, false) > 0;
	bool have_upper = upper && p(*upper, parameters, false, false) > 0;

	if (have_lower && have_upper) {
	    if (p(ll, parameters, false, false) < 0.5) {
		//Use upper tail
		loglik -= log(p(ll, parameters, false, false) -
			      p(*upper, parameters, false, false));
	    }
	    else {
		//Use lower tail
		loglik -= log(p(*upper, parameters, true, false) - 
			      p(ll, parameters, true, false));
	    }
	}
	else if (have_lower) {
	    loglik -= p(ll, parameters, false, true);
	}
	else if (have_upper) {
	    loglik -= p(*upper, parameters, true, true);
	}
    }

    return loglik;
}


double 
RScalarDist::randomSample(vector<double const *> const &parameters,
			  double const *lower, double const *upper,
			  RNG *rng) const
{
    if (lower || upper) {

	double plower = 0, pupper = 1;
	if (lower) {
	    plower = calPlower(*lower, parameters);
	}
	if (upper) {
	    pupper = calPupper(*upper, parameters);
	}
	
	double u = plower + rng->uniform() * (pupper - plower);
	return q(u, parameters, true, false);
    }
    else {
	return r(parameters, rng);
    }

}

bool RScalarDist::canBound() const
{
    return true;
}

bool RScalarDist::isDiscreteValued(vector<bool> const &mask) const
{
    return _discrete;
}

bool RScalarDist::discrete() const
{
    return _discrete;
}

unsigned int RScalarDist::npar() const
{
    return _npar;
}
