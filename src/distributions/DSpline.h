/*
##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, National Council On Compensation Insurance Inc.,                        ##
##                                                                                              ##
##    This file is part of lossDev.                                                             ##
##                                                                                              ##
##    lossDev is free software: you can redistribute it and/or modify                           ##
##    it under the terms of the GNU General Public License as published by                      ##
##    the Free Software Foundation, either version 3 of the License, or                         ##
##    (at your option) any later version.                                                       ##
##                                                                                              ##
##    This program is distributed in the hope that it will be useful,                           ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of                            ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             ##
##    GNU General Public License for more details.                                              ##
##                                                                                              ##
##    You should have received a copy of the GNU General Public License                         ##
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.                     ##
##                                                                                              ##
##################################################################################################
*/

#ifndef DSPLINE_H_
#define DSPLINE_H_

#include <distribution/ArrayDist.h>

class DSpline : public ArrayDist
{
public:
	DSpline();
	virtual ~DSpline();
	    /**
	     * @param x Value at which to evaluate the likelihood.
	     *
	     * @param parameters Vector of parameter values of the
	     * distribution.
	     * 
	     * @param dims Vector of dimensions of the parameters.
	     *
	     * @param lbound Lower bound for truncated distributions. If the
	     * distribution is not truncated then this should be a NULL pointer.
	     *
	     * @param ubound Upper bound for truncated distributions. If the
	     * distribution is not truncated then this should be a NULL pointer.
	     * 
	     * @returns the log likelihood.  If the likelihood should be zero
	     * because x is inconsistent with the parameters then -Inf is
	     * returned. If the parameters are invalid
	     * (i.e. checkParameterValue returns false), then the return value
	     * is undefined.
	     *
             * NOTE: This is currently only valid for the 7th and 8th parameters
	     */
	    virtual double 
		logLikelihood(double const *x, unsigned int length,
			      std::vector<double const *> const &parameters,
			      std::vector<std::vector<unsigned int> > const &dims,
			      double const *lbound, double const *ubound)
		const;
	    /**
	     * Draws a random sample from the distribution. 
	     *
	     * @param x Array to which the sample values are written
	     *
	     * @param parameters  Vector of parameter values at which
	     * to evaluate the likelihood. This vector should be of length
	     * npar().
	     *
	     * @param lbound Lower bound, for truncated distributions, or a NULL
	     * pointer if the distribution is not truncated.
	     *
	     * @param ubound Upper bound, for truncated distributions, or a NULL
	     * pointer if the distribution is not truncated.

	     * @param rng pseudo-random number generator to use.
	     *
	     * @exception length_error 
	     */
	    virtual void 
		randomSample(double *x, unsigned int length,
			     std::vector<double const *> const &parameters,
			     std::vector<std::vector<unsigned int> > const  &dims,
			     double const *lbound, double const *ubound, RNG *rng) 
		const;
	    /**
	     * Returns a typical value from the distribution.  The meaning of
	     * this will depend on the distribution, but it will normally be a
	     * mean, median or mode.
	     *
	     * @param x Array to which the sample values are written
	     *
	     * @param parameters  Vector of parameter values at which
	     * to evaluate the likelihood. This vector should be of length
	     * npar().
	     *
	     * @param dims Vector of parameter dimensions.
	     *
	     * @param lbound Lower bound, for truncated distributions, or a NULL
	     * pointer if the distribution is not truncated.
	     *
	     * @param ubound Upper bound, for truncated distributions, or a NULL
	     * pointer if the distribution is not truncated.

	     *
	     * @exception length_error 
	     */
	    virtual void 
		typicalValue(double *x, unsigned int length,
			     std::vector<double const *> const &parameters,
			     std::vector<std::vector<unsigned int> > const &dims,
			     double const *lbound, double const *ubound)
		const;

	    /**
	     * The lower and upper limits on the support of the distribution,
	     * given the parameters. 
	     *
	     * @param lower Pointer to the start of an array to which the lower
	     * limit will be written.
	     *  
	     * @param upper Pointer to the start of an array to which the lower
	     * limit will be written.
	     *
	     * @param length Length of the lower and upper arrays
	     *
	     * @param parameters Vector of parameters at which to evaluate the
	     * support
	     *
	     * @param dims vector of dimensions of the parameters
	     */
	    virtual void support(double *lower, double *upper, unsigned int length,
				 std::vector<double const *> const &parameters,
				 std::vector<std::vector<unsigned int> > const &dims) 
		const;
	    /**
	     * Indicates whether the support of the distribution is fixed.
	     *
	     * @param fixmask Boolean vector of length npar() indicating which
	     * parameters have fixed values.
	     */
	    virtual bool isSupportFixed(std::vector<bool> const &fixmask) const;
	    /**
	     * Checks that dimensions of the parameters are correct.
	     *
	     * This function only needs to be run once for each parameter
	     * vector. Thereafter, the values of the parameters will change,
	     * but the dimensions will not.
	     */
	    virtual bool 
		checkParameterDim (std::vector<std::vector<unsigned int> > const &parameters) 
		const;
	    /**
	     * Checks that the values of the parameters are consistent with
	     * the distribution. For example, some distributions require
	     * than certain parameters are positive, or lie in a given
	     * range.
	     *
	     * This function assumes that checkParameterDim returns true.
	     */
	    virtual bool 
		checkParameterValue(std::vector<double const *> const &parameters,
				    std::vector<std::vector<unsigned int> > const &dims) const;
	    /**
	     * Calculates what the dimension of the distribution should be,
	     * based on the dimensions of its parameters. Since the majority
	     * of distributions are scalar valued, a default function is provided
	     * that returns a vector of length 1 with element 1. This must be
	     * overloaded by multivariate distributions.
	     */
	    virtual std::vector<unsigned int> 
		dim (std::vector <std::vector<unsigned int> > const &args) const;
	    /**
	     * Returns the number of degrees of freedom of the distribution
	     * given the dimensions of the parameters. By default this is the
	     * product of the elements of the dimension vector returned by
	     * Distribution#dim. However, some distributions are constrained:
	     * and the support occupies a lower dimensional subspace. In this
	     * case, the df member function must be overrideen.
	     */
	    //virtual unsigned int df(std::vector<std::vector<unsigned int> > const &dims) const;
	    /**
	     * Tests for a location parameter.  A parameter of a distribution
	     * is considered to be a location parameter if, when it's value is
	     * incremented by X, the whole distribution is shifted by X,
	     * indpendently of the other parameter values.
	     * 
	     * This is a virtual function, for which the default implementation
	     * always returns false. Distributions with location parameters must
	     * overload this function.
	     *
	     * @param index Index number (starting from 0) of the parameter to be
	     * tested.
	     */
	    //virtual bool isLocationParameter(unsigned int index) const;
	    /**
	     * Tests for a scale parameter.  A parameter of a distribution is
	     * considered to be a scale parameter if, when it's value is
	     * multiplied by X, the whole distribution multiplied by X,
	     * indpendently of the other parameter values.
	     * 
	     * Note that this definition excludes "location-scale" models:
	     * i.e. if the density of y takes the form (1/b)*f((y-a)/b) then b
	     * is not considered a scale parameter.
	     *
	     * This is a virtual function, for which the default
	     * implementation always returns false. Distributions with scale
	     * parameters must overload this function.
	     *
	     * @param index Index number (starting from 0) of the parameter to
	     * be tested.
	     */
	    //virtual bool isScaleParameter(unsigned int index) const;
};

#endif /*DSPLINE_H_*/
