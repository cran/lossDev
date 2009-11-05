#ifndef SLICER_H_
#define SLICER_H_

#include <sampler/DensityMethod.h>

class Sampler;

/**
 * @short Slice Sampling updating method
 * 
 * The slice sampler is defined by Neal R (2003) Slice Sampling,
 * Ann. Statist.  31; 3: 705-767.
 *
 * The Slicer class is written with a high level of abstraction to
 * allow it to be easily adapted: for example, a slice sampler can be
 * used to sample from a discrete distribution or from a univariate
 * "slice" of a multivariate distribution.
 */
class Slicer : public DensityMethod 
{
    double _width;
    bool _adapt;
    unsigned int _max;
    double _sumdiff;
    unsigned int _iter;
    bool accept(double xold, double xnew, double z, double L, double R,
		double lower, double upper);
public:
    /**
     * Constructor for Slice Sampler
     *
     * @param width Initial width of slice
     *
     * @param max Maximum number of times initial width of slice will
     * increase at each iteration.
     */
    Slicer(double width, unsigned int max);
    /**
     * Update the current value using the "stepping" method. A Sampler
     * that uses the Slicer DensityMethod can implement Sampler#update
     * by calling this function.
     */
    void updateStep(RNG *rng);
    /**
     * Update the current value using the "doubling" method. A Sampler
     * that uses the Slicer DensityMethod can implement Sampler#update
     * by calling this function.
     */
    void updateDouble(RNG *rng);
    /**
     * Returns the current value of the Slicer. A Slicer may have it's
     * own value, distinct from the value of the StochasticNode being
     * updated.
     */
    virtual double value() const = 0;
    /**
     * Sets the value of the sampler, and additionally updates all the
     * deterministic descendants. This function is virtual because the
     * value of the sampler and the value of the sampled node must be
     * defined by a subclass of Slicer.
     *
     * @param x value to set
     */
    virtual void setValue(double x) = 0;
    /**
     * Gets the lowest and highest possible values of the sampler
     */
    virtual void getLimits(double *lower, double *upper) const = 0;
    /**
     * Turns off adaptive mode.  The current adaptation test is very
     * basic, and will return true if a minimum number of iterations
     * (50) have taken place.
     */
    bool adaptOff();
    /**
     * The slicer method is adaptive.  The step size adapts to the
     * mean distance between consecutive updates
     */
   bool isAdaptive() const;
};

#endif /* SLICER_H_ */
