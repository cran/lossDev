#ifndef DSUM_H_
#define DSUM_H_

#include <distribution/DistScalar.h>

/**
 * @short Sum of two discrete random variables
 */
class DSum : public DistScalar {
public:
    DSum();

    double scalarLogLikelihood(double x, 
			       std::vector<double const *> const &parameters,
			       double const *lower, double const *upper) const;
    double scalarRandomSample(std::vector<double const *> const &parameters,
			      double const *lower, double const *upper,
			      RNG *rng) const;
    bool checkParameterDiscrete (std::vector<bool> const &mask) const;
    double l(std::vector<double const *> const &parameters) const;
    double u(std::vector<double const *> const &parameters) const;
    double typicalScalar(std::vector<double const *> const &parameters,
			 double const *lower, double const *upper) const;
    bool isSupportFixed(std::vector<bool> const &fixmask) const;
    unsigned int df(std::vector<std::vector<unsigned int> > const &dims) const;
};

#endif /* DSUM_H_ */
