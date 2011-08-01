#ifndef DUNIFOV_H_
#define DUNIFOV_H_

#include <distribution/ScalarDist.h>

/**
 * <pre>
 * x ~ dunif(a, b)
 * f(x|a,b) = 1/(a - b) ; a <= x <= b
 * </pre>
 * @short Uniform distribution
 */
class DUnifOV : public ScalarDist {
  public:
    DUnifOV();
    double logDensity(double x, PDFType type,
                      std::vector<double const *> const &parameters,
                      double const *lower, double const *upper) const;
    double randomSample(std::vector<double const *> const &parameters,
			double const *lower, double const *upper,
			RNG *rng) const;
    double typicalValue(std::vector<double const *> const &parameters,
			double const *lower, double const *upper) const;
    double l(std::vector<double const*> const &parameters) const;
    double u(std::vector<double const*> const &parameters) const;
    /** 
     * Checks that a < b
     */
    bool checkParameterValue(std::vector<double const *> const &parameters)
	const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

#endif /* DUNIFOV_H_ */
