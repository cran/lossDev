#ifndef DT_OV_H_
#define DT_OV_H_

#include <distribution/DistScalarRmath.h>

/**
 * t-distribution on k degrees of freedom, with median mu and
 * scale parameter tau.
 * <pre>
 * f(x|mu, tau, k)
 * f(x|0,1,k) = Gamma((k+1)/2) / (sqrt(k*pi) Gamma(k/2)) (1 + x^2/k)^-((k+1)/2)
 * </pre>
 * @short t distribution
 */
class DTOV : public DistScalarRmath {
 public:
  DTOV();

  double d(double x, std::vector<double const *> const &parameters, 
	   bool log) const;
  double p(double x, std::vector<double const *> const &parameters, bool lower,
	   bool log) const;
  double q(double x, std::vector<double const *> const &parameters, bool lower,
	   bool log) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Check that tau > 0 and k > 0
   */
  bool checkParameterValue (std::vector<double const *> const &parameters,
			    std::vector<std::vector<unsigned int> > const &dims)
    const;
};

#endif /* DT_OV_H_ */
