#ifndef DDIRCH_H_
#define DDIRCH_H_

#include <distribution/Distribution.h>

struct RNG;

/**
 * @short Dirichlet distribution
 *
 * Zero shape parameters are allowed.  These represent structural
 * zeros: when x ~ ddirch(alpha) is forward sampled, x[i] = 0 when
 * alpha[i] = 0. To avoid trapping states in the model, structural
 * zeros are only allowed when the array of shape parameters is
 * fixed.
 *
 * <pre>
 * p[] ~ ddirch(alpha[])
 * f(p | alpha) = C * prod(p^alpha)
 * </pre>
 */
class DDirch : public Distribution {
public:
    DDirch();

    double logLikelihood(double const *x, unsigned int length,
			 std::vector<double const *> const &parameters,
			 std::vector<std::vector<unsigned int> > const &dims,
			 double const *lower, double const *upper) const;
    void randomSample(double *x, unsigned int length,
		      std::vector<double const *> const &parameters,
		      std::vector<std::vector<unsigned int> > const &dims,
		      double const *lower, double const *upper, RNG *rng) const;
    void typicalValue(double *x, unsigned int length,
		      std::vector<double const *> const &parameter,
		      std::vector<std::vector<unsigned int> > const &dims,
		      double const *lower, double const *upper) const;
    std::vector<unsigned int> 
	dim(std::vector<std::vector<unsigned int> > const &dims) const;
    /**
     * Checks that alpha is a vector of length at least 2
     */
    bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) const;
    /**
     * Checks that each element of alpha is >= 0.
     *
     * Structural zeros are allowed in the Dirichlet distribution.
     * These are represented by the elements of alpha that are set to
     * zero.  This is permitted only if alpha is fixed and there is at
     * least one non-zero element of alpha.
     */
    bool checkParameterValue(std::vector<double const *> const &parameters,
			     std::vector<std::vector<unsigned int> > const &dims)
	const;
    void support(double *lower, double *upper, unsigned int length,
		 std::vector<double const *> const &parameters,
		 std::vector<std::vector<unsigned int> > const &dims) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    unsigned int df(std::vector<std::vector<unsigned int> > const &dims) const;
};

#endif /* DDIRCH_H_ */
