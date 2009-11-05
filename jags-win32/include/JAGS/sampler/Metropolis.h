#ifndef METROPOLIS_H_
#define METROPOLIS_H_

#include <sampler/DensityMethod.h>
#include <vector>

class StochasticNode;

/**
 * @short Metropolis-Hastings update method
 *
 * This class is used by Metropolis Hastings samplers.  It provides
 * only basic infrastructure.
 *
 * The Metropolis Hasting update method has its own value, which is
 * distinct from the value of the sampled nodes.  All member functions
 * of the Metropolis class work in terms of this value.  This is done
 * to allow parameter transformations: for example, if a node only
 * takes positive values, it may be more efficient for the sampler to
 * work on log-transformed values of the sampled node. Parameter
 * transformations are defined by the member functions
 * Metropolis#transform and Metropolis#untransform.
 *
 * A Subclass of Metropolis must provide an implementation of 
 * the virtual functions Metropolis#setValue and Metropolis#rescale.
 *
 * The Metropolis class provides no update member function. A subclass
 * of Metropolis must provide this. It should contain one or more
 * calls to Metropolis#propose, then calculate the acceptance
 * probability, and then call the function Metropolis#accept.
 */
class Metropolis : public DensityMethod
{
    bool _adapt;
    double *_value;
    double *_last_value;
    unsigned int _length;
    Metropolis(Metropolis const &);
    Metropolis &operator=(Metropolis const &);
public:
    Metropolis(std::vector<StochasticNode *> const &nodes);
    ~Metropolis();
    /**
     * Sets the initial value of the Metropolis object by taking
     * the current value of the sampled nodes and transforming them with
     * the Metropolis#untransform function.
     *
     * Initialization cannot be done when the Metropolis object is
     * constructed, as it depends on the virtual untransform function.
     *
     * A subclass of Metropolis must not overload this function. It
     * should overload initMetropolis instead. 
     */
    void initialize(DensitySampler *sampler, unsigned int chain);
    /**
     * Extra initialization function called by Metropolis#initialize.
     * A subclass of Metropolis should overload this function instead of
     * the initialize function if it needs additional initialization of
     * private data. The default function does nothing.
     */
    virtual void initMetropolis(DensitySampler *sampler, unsigned int chain);
    /**
     * Returns the current value array of the Metropolis object. 
     */
    double const *value() const;
    /**
     * Returns the length of the value array. 
     */
    unsigned int value_length() const;
    /**
     * Sets the value of the Metropolis object. 
     *
     * @param value Pointer to the beginning of an array of values
     *
     * @param length Length of the supplied value array
     */
    void propose(double const *value, unsigned int length);
    /**
     * Accept current value with probabilty p. If the current value is
     * not accepted, the Metropolis object reverts to the value at the
     * last successful call to accept.
     *
     * @param rng Random number generator.
     *  
     * @param p Probability of accepting the current value.
     *
     * @returns success indicator
     */
    bool accept(RNG *rng, double p);
    /**
     * Rescales the proposal distribution. This function is called by
     * Metropolis#accept when the sampler is in adaptive
     * mode. Rescaling may depend on the acceptance
     * probability.
     *
     * @param p Acceptance probability
     */
    virtual void rescale(double p) = 0;
    /**
     * The Metropolis-Hastings method is adaptive. The process of
     * adaptation is specific to each subclass and is defined by the
     * rescale member function
     */
    bool isAdaptive() const;
    /**
     * Turns off adaptive mode
     */
    bool adaptOff();
    /**
     * Tests whether adaptive mode has been successful (e.g. by testing
     * that the acceptance rate lies in an interval around the target 
     * value). This function is called by Metropolis#adaptOff;
     */
    virtual bool checkAdaptation() const = 0;
    /**
     * Transforms the value of the Metropolis object to an array of
     * concatenated values for the sampled nodes.
     *
     * This function is called by Metropolis#propose
     */
    virtual void transform(double const *v, unsigned int length,
			   double *nv, unsigned int nlength) const = 0;
    /**
     * Transforms the concatenated values of the sampled nodes to 
     * a value for the Metropolis object.
     *
     * This function is called by Metropolis#initialize
     */
    virtual void untransform(double const *nv, unsigned int nlength,
			     double *v, unsigned int length) const = 0;
};

#endif /* METROPOLIS_H_ */
