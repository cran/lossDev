#ifndef DBERN_H_
#define DBERN_H_

#include <distribution/DistScalarRmath.h>

/**
 * @short Bernoulli distribution
 * <pre>
 * R ~ dbern(p)
 * f(r | p) = p^r * (1 - p)^(1 -r) ; r in 0:1
 * </pre>
 */
class DBern : public DistScalarRmath {
public:
    DBern();

    double d(double x, std::vector<double const *> const &parameters, 
	     bool give_log) const;
    double p(double x, std::vector<double const *> const &parameters, 
	     bool lower,
	     bool give_log) const;
    double q(double p, std::vector<double const *> const &parameters, 
	     bool lower,
	     bool log_p) const;
    double r(std::vector<double const *> const &parameters, RNG *rng) const;
    /** Checks that p lies in the open interval (0,1) */
    bool 
      checkParameterValue(std::vector<double const *> const &parameters,
			  std::vector<std::vector<unsigned int> > const &dims)
      const;
};

#endif /* DBERN_H_ */
