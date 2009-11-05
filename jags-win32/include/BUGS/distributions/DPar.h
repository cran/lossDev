#ifndef DPAR_H_
#define DPAR_H_

#include <distribution/DistScalarRmath.h>

/**
 * <pre>
 * x ~ dpar(alpha, c);
 * f(x|alpha,c) = alpha * c^alpha * x^-(alpha + 1); x > c
 * </pre>
 * @short Pareto distribution
 */
class DPar : public DistScalarRmath {
 public:
  DPar();

  double d(double x, std::vector<double const *> const &parameters, bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  double l(std::vector<double const*> const &parameters) const;
  double u(std::vector<double const*> const &parameters) const;
  /** 
   * Checks that alpha > 0, c > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned int> > const &dims) 
    const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

#endif /* DPAR_H_ */
