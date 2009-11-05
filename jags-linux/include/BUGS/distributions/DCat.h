#ifndef DCAT_H_
#define DCAT_H_

#include <distribution/Distribution.h>

/**
 * @short Categorical distribution
 * <pre>
 * R ~ dcat(p[])
 * f(r|p[]) = p[r] ; r in 1:dim(p)
 * </pre>
 */
class DCat : public Distribution {
public:
    DCat();

    double logLikelihood(double const *x, unsigned int length,
			 std::vector<double const *> const &parameters,
			 std::vector<std::vector<unsigned int> > const &dims,
			 double const *lower, double const *upper) const;
    void randomSample(double *x, unsigned int length,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned int> > const &dims,
		      double const *lbound, double const *ubound,
		      RNG *rng) const;
    void typicalValue(double *x, unsigned int length,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned int> > const &dims,
		      double const *lbound, double const *ubound) const;
    /**
     * Checks that  p is a vector of length at least 2
     */
    bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) 
	const;
    /**
     * Checks that all elements of p are positive
     */
    bool 
	checkParameterValue(std::vector<double const*> const &parameters,
			    std::vector<std::vector<unsigned int> > const &dims)
	const;
    void support(double *lower, double *upper, unsigned int length,
		 std::vector<double const *> const &parameters,
		 std::vector<std::vector<unsigned int> > const &dims) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

#endif /* DCAT_H_ */
