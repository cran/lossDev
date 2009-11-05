#ifndef DLNORM_H_
#define DLNORM_H_

#include <distribution/DistScalarRmath.h>

/** 
 * @short Lognormal distribution
 * <pre>
 * x ~ dlnorm(mu, tau)
 * f(x|mu,tau) = (sqrt(tau)/x) * exp(tau * (log(x) - mu)^2 / 2)
 * </pre>
 */
class DLnorm : public DistScalarRmath {
 public:
  DLnorm();

  double d(double x, std::vector<double const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Checks that tau > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned int> > const &dims)
    const;
};

#endif /* DLNORM_H_ */
