#ifndef DGAMMA_H_
#define DGAMMA_H_

#include <distribution/DistScalarRmath.h>

/**
 * @short gamma distribution
 * <pre>
 * X ~ dgamma(r, mu)
 * f(x|r,mu) = mu^r * x^(r - 1) * exp(-mu * x)
 * </pre>
 */
class DGamma : public DistScalarRmath {
 public:
  DGamma();

  double d(double x, std::vector<double const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Checks that r > 0, mu > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned int> > const & dims)
    const;
  /**
   * The mean
   */
  double  typicalScalar(std::vector<double const *> const &par,
			double const *lower, double const *upper) const;
};

#endif /* DGAMMA_H_ */
