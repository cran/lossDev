#ifndef DMULTI_H_
#define DMULTI_H_

#include <distribution/Distribution.h>

struct RNG;

/**
 * <pre>
 * X[] ~ dmulti(p[], N)
 * f(x | p, N) = prod (p^x) ; sum(x) = N
 * </pre>
 * @short Multinomial distribution
 */
class DMulti : public Distribution {
public:
  DMulti();

  double logLikelihood(double const *x, unsigned int length, 
		       std::vector<double const *> const &parameters,
                       std::vector<std::vector<unsigned int> > const &dims,
		       double const *lower, double const *upper) const;
  void randomSample(double *x, unsigned int length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const &dims,
		    double const *lower, double const *upper, RNG *rng) const;
  void typicalValue(double *x, unsigned int length,
		    std::vector<double const *> const &par,
		    std::vector<std::vector<unsigned int> > const &dims,
		    double const *lower, double const *upper) const;
  /**
   * Checks that elements of p lie in range (0,1) and 
   * and sum to 1. Checks that N >= 1
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
                           std::vector<std::vector<unsigned int> > const &dims)
      const;
  /** Checks that p is a vector and N is a scalar */
  bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dim) const;
  /** Checks that N is discrete-valued */
  bool checkParameterDiscrete(std::vector<bool> const &mask) const;
  std::vector<unsigned int> 
      dim(std::vector<std::vector<unsigned int> > const &dim) const;
  void support(double *lower, double *upper, unsigned int length,
	       std::vector<double const *> const &parameters,
	       std::vector<std::vector<unsigned int> > const &dims) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
  unsigned int df(std::vector<std::vector<unsigned int> > const &dims) const;
};

#endif /* DMULTI_H_ */
