#ifndef DWISH_H_
#define DWISH_H_

#include <distribution/Distribution.h>

/**
 * <pre>
 * x[] ~ dwish(R[,], k)
 * </pre>
 * @short Wishart distribution
 */
class DWish : public Distribution {
public:
  DWish();

  double logLikelihood(double const *x, unsigned int length,
		       std::vector<double const *> const &parameters,
		       std::vector<std::vector<unsigned int> > const &dims,
		       double const *lower, double const *upper) const;
  void randomSample(double *x, unsigned int length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const &dims,
		    double const *lower, double const *upper, RNG *rng) const;
  //FIXME: Can we retire this?
  static void randomSample(double *x, int length,
                           double const *R, double k, int nrow,
                           RNG *rng);
  void typicalValue(double *x, unsigned int length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const &dims,
		    double const *lower, double const *upper) const;
  /**
   * Checks that R is a square matrix and k is a scalar
   */
  bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) 
      const;
  /**
   * Checks that R is symmetric and k >= nrow(R). There is
   * currently no check that R is positive definite
   */
  bool checkParameterValue(std::vector<double const *> const &parameters,
			   std::vector<std::vector<unsigned int> > const &dims)
      const;
  std::vector<unsigned int> 
      dim(std::vector<std::vector<unsigned int> > const &dims) const;
  void support(double *lower, double *upper, unsigned int length,
	       std::vector<double const *> const &parameters,
	       std::vector<std::vector<unsigned int> > const &dims) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
  unsigned int df(std::vector<std::vector<unsigned int> > const &dims) const;
};

#endif /* DWISH_H_ */
