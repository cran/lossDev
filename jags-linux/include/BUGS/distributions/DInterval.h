#ifndef DINTERVAL_H_
#define DINTERVAL_H_

#include <distribution/Distribution.h>

/**
 * @short Interval censored distribution
 * <pre>
 * i ~ dinterval(t, cutpoints[])
 * f(i|t) = 1 if t < cutpoints[i] and t >= cutpoints[i-1]
 *        = 0 otherwise
 * </pre>
 */
class DInterval : public Distribution {
public:
    DInterval();
  
    double logLikelihood(double const *x, unsigned int length,
			 std::vector<double const *> const &parameters,
			 std::vector<std::vector<unsigned int> > const &dims,
			 double const *lower, double const *upper) const;

    void randomSample(double *x, unsigned int length,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned int> > const &dims,
		      double const *lower, double const *upper, RNG *rng) const;
    void typicalValue(double *x, unsigned int length,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned int> > const &dims,
		      double const *lower, double const *upper) const;

    /**
     * Checks that t is scalar and that cutpoints is a vector
     */
    bool checkParameterDim(std::vector<std::vector<unsigned int> > const &parameters) const;
    /**
     * Checks that cutpoints are in ascending order
     */
    bool checkParameterValue(std::vector<double const *> const &parameters,
			     std::vector<std::vector<unsigned int> > const &dims) const;
    void support(double *lower, double *upper, unsigned int length,
		 std::vector<double const *> const &parameters,
		 std::vector<std::vector<unsigned int> > const &dims) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    unsigned int df(std::vector<std::vector<unsigned int> > const &dims) const;

};

#endif /* DINTERVAL_H_ */
