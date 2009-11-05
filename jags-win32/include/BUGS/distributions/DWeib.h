#ifndef DWEIB_H_
#define DWEIB_H_

#include <distribution/DistScalarRmath.h>

/**
 * <pre>
 * x ~ dweib(a, b)
 * f(x|b,a) = a * b * x^(b - 1) * exp (- a * x^b)
 * </pre>
 * @short Weibull distribution
 */
class DWeib : public DistScalarRmath {
public:
  DWeib();

  double d(double x, std::vector<double const *> const &parameters, 
	   bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /** 
   * Checks that a > 0, b > 0
   */
  bool checkParameterValue (std::vector<double const *> const &parameters,
			    std::vector<std::vector<unsigned int> > const &dims)
    const;

};

#endif /* DWEIB_H_ */
