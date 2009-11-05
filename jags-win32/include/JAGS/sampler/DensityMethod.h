#ifndef DENSITY_METHOD_H_
#define DENSITY_METHOD_H_

#include <string>

class DensitySampler;
class RNG;

/**
 * @short Updates a Sampler using the log density of the sampled nodes.
 *
 * Updates a single chain of a DensitySampler. A typical
 * DensityMethod uses only the log density of the sampled nodes by
 * calling Sampler#logFullConditional.
 */
class DensityMethod
{
protected:
    DensitySampler *_sampler;
    unsigned int _chain;
public:
    DensityMethod();
    virtual ~DensityMethod();
    /**
     * Sets the values of the protected data members _sampler and _chain
     * so that subclasses of DensityMethod can access them. This
     * function is called by the constructor of DensitySampler.
     */
    void setData(DensitySampler *sampler, unsigned int chain);
    /**
     * A subclass of DensityMethod may optionally set data members
     * when this function is called. The default does nothing.
     */
    virtual void initialize(DensitySampler *sampler,  unsigned int chain);
    /**
     * Updates the sampler and chain defined by the last call to setData
     */
    virtual void update(RNG *rng) = 0;
    /**
     * Turns off adaptive mode, returning true if an adaptation test is
     * passed.
     *
     * @see Sampler#adaptOff
     */
    virtual bool adaptOff() = 0;
    /**
     * Indicates whether the update method has an adaptive mode
     *
     * @see Sampler#isAsaptive
     */
    virtual bool isAdaptive() const = 0;
    /**
     * Returns an informative name for the update method
     *
     * @see Sampler#name
     */
    virtual std::string name() const = 0;
};
    
#endif /* DENSITY_METHOD_H_ */
