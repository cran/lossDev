#ifndef DPOIS_H_
#define DPOIS_H_

#include <distribution/DistScalarRmath.h>

/**
 * <pre>
 * r ~ dpois(lambda)
 * f(r|lambda) = exp(-lambda) * lambda^r / r! ; r = 0, 1 ,...
 * </pre>
 * @short Poisson distribution
 */
class DPois : public DistScalarRmath {
 public:
  DPois();

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

#endif /* DPOIS_H_ */
