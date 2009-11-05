#ifndef DNEGBIN_H_
#define DNEGBIN_H_

#include <distribution/DistScalarRmath.h>

/**
 * <pre>
 * x ~ dnegbin(p, r)
 * f(x|p,r) = ((x+r-1)!/(x!*(r-1)!)) * p^r * (1-p)^x
 * </pre>
 * @short Negative Binomial distribution
 */
class DNegBin : public DistScalarRmath {
 public:
  DNegBin();

  double d(double x, std::vector<double const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /** Checks that r is discrete-valued */
  bool checkParameterDiscrete (std::vector<bool> const &mask) const;
  /**
   * Checks that p lies in the interval (0,1) and r > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned int> > const &dims)
    const;

};

#endif /* DNEGBIN_H_ */
