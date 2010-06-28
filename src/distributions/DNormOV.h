#ifndef DNORMOV_H_
#define DNORMOV_H_

#include "RScalarDist.h"

/**
 * <pre>
 * x ~ dnorm(mu, tau)
 * f(x | mu, tau) = sqrt(tau) * exp(-1/2 * tau * (x - mu)^2)
 * </pre>
 * @short Normal distribution
 */
class DNormOV : public RScalarDist {
 public:
  DNormOV();

  double d(double x, std::vector<double const *> const &parameters, 
	   bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Checks that tau > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters) const;
  /**
   * Exploits the capacity to sample truncted normal distributions
   * that is built into the JAGS library, overloading the generic
   * functionality of RScalarDist.
   */
  double randomSample(std::vector<double const *> const &par,
		      double const *lower, double const *upper,
		      RNG *rng) const;
};

#endif /* DNORMOV_H_ */
