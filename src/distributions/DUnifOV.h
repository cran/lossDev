#ifndef DUNIF_OV_H_
#define DUNIF_OV_H_

#include <distribution/DistScalarRmath.h>

/**
 * <pre>
 * x ~ dunif(a, b)
 * f(x|a,b) = 1/(a - b) ; a <= x <= b
 * </pre>
 * @short Uniform distribution
 */
class DUnifOV : public DistScalarRmath {
 public:
  DUnifOV();

  double d(double x, std::vector<double const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  double l(std::vector<double const*> const &parameters) const;
  double u(std::vector<double const*> const &parameters) const;
  /** 
   * Checks that a < b
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned int> > const &dims) 
    const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

#endif /* DUNIF_OV_H_ */
