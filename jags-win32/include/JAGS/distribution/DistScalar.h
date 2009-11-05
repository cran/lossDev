#ifndef DIST_SCALAR_H_
#define DIST_SCALAR_H_

#include <distribution/Distribution.h>

struct RNG;

/**
 * Enumerates three possible ranges of support for a scalar random
 * variable.
 *
 * DIST_UNBOUNDED for support on the whole real line
 *
 * DIST_POSITIVE for support on values > 0
 *
 * DIST_PROPORTION for support on values in [0,1]
 *
 * DIST_SPECIAL for other distributions, e.g. distributions where the
 * support depends on the parameter values.
 */
enum Support {DIST_UNBOUNDED, DIST_POSITIVE, DIST_PROPORTION, DIST_SPECIAL};

/**
 * Base class for scalar valued distributions, whose parameters are
 * also scalars.
 * *
 * @short Scalar distributions
 */
class DistScalar : public Distribution
{
  const Support _support;
 public:
  /**
   * Constructor
   *
   * @param name BUGS language name of distribution
   *
   * @param npar Number of parameters, excluding upper and lower bound
   *
   * @param support Support of distribution
   *
   * @param boundable Logical flag indicating whether the distribution
   * can be bounded using the T(,) construct.
   *
   * @param discrete Logical flag indicating whether the distribution
   * is discrete-valued.
   */
  DistScalar(std::string const &name, unsigned int npar,
	     Support support, bool canbound, bool discrete);
  double logLikelihood(double const *x, unsigned int length,
		       std::vector<double const *> const &parameters,
		       std::vector<std::vector<unsigned int> > const  &dims,
                       double const *lower, double const *upper)
      const;
  void randomSample(double *x, unsigned int length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const  &dims,
                    double const *lower, double const *upper,
		    RNG *r) const;
  void typicalValue(double *x, unsigned int length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const &dims,
                    double const *lower, double const *upper) const;
  /**
   * Checks that parameters are scalar
   */
  bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims)
     const;
  /**
   * This implementation of calculates the support
   * based on the DistScalar#l and DistScalar#u functions.
   */
  void support(double *lower, double *upper, unsigned int length, 
	       std::vector<double const *> const &parameters,
	       std::vector<std::vector<unsigned int> > const &dims) const;
  /**
   * This implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
  /**
   * Lower limit of distribution, given parameters.  If the
   * distribution has no lower limit, this should return JAGS_NEGINF.
   *
   * The default implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  virtual double l(std::vector<double const *> const &parameters) const;
  /**
   * Upper limit of distribution, given parameters. If the
   * distribution has no upper limit, this should return JAGS_POSINF.
   *
   * The default implementation should be used for distributions with
   * Support DIST_UNBOUNDED, DIST_POSITIVE and DIST_PROPORTION. If
   * the Support is DIST_SPECIAL, this must be overloaded.
   */
  virtual double u(std::vector<double const *> const &parameters) const;
  /**
   * Simplified version of loglikelihood function for scalar distributions
   */
  virtual double 
      scalarLogLikelihood(double x, 
			  std::vector<double const *> const &parameters,
			  double const *lower, double const *upper)
      const = 0;
  /**
   * Simplified version of the randomSample function for scalar distributions
   */
  virtual double 
    scalarRandomSample(std::vector<double const *> const &parameters, 
		       double const *lower, double const *upper,
		       RNG *rng) const = 0;
  /**
   * Simplified versino of typicalValue function for scalar distributions
   */
  virtual double typicalScalar(std::vector<double const *> const &parameters,
			       double const *lower, double const *upper)
    const = 0;
  /**
   * Simplified version of support function for scalar distributions
   */
  void support(double *lower, double *upper, 
	       std::vector<double const *> const &parameters) const;
};

#endif /* DIST_SCALAR_H_ */
