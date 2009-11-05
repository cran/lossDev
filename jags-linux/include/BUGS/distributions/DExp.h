#ifndef DEXP_H_
#define DEXP_H_

#include <distribution/DistScalarRmath.h>

/**
 * @short Exponential distribution
 * <pre>
 *  x ~ dexp(lambda)
 *  f(x | lambda) = lambda * exp(-lambda * x) ; x >= 0
 * </pre>
 */
class DExp : public DistScalarRmath {
 public:
  DExp();

  double d(double x, std::vector<double const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Checks that lambda > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned int> > const &dims)
    const;
};

#endif /* DEXP_H_ */
