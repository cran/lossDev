#ifndef REAL_SLICER_OV_H_
#define REAL_SLICER_OV_H_

#include <sampler/SampleMethod.h>
#include <sampler/GraphView.h>

class StochasticNode;


/**
 * Slice sampler for real-valued distributions
 */
class RealSlicerOV : public SampleMethod
{
    GraphView const *_gv;
    unsigned int _chain;

    double _probOfOverrelaxed;
    unsigned int _endpointAccuracy;

    unsigned int _maxOV;
    double _widthOV;
    unsigned int _iterOV;
    double _sumdiffOV;

    bool _adaptOV;
 


public:
    /**
     * Constructor for Slice Sampler
     * @param node Node to sample
     * @param width Initial width of slice
     * @param maxwidth Maximal width of slice as a multiple of the width
     * parameter
     * @param nburn Length of burnin
     */
    RealSlicerOV(GraphView const * gv, unsigned int chain);
    double value() const;
    void setValue(double value) const;
    void getLimits(double *lower, double *upper) const;
    virtual void update(RNG *rng);
    void updateOverrelaxed(RNG *rng);
    void updateStep(RNG *rng);
    virtual std::string name() const;
    static bool canSample(StochasticNode const *node);
    virtual bool isAdaptive() const;
    virtual bool adaptOff();
    double logFullConditional() const;
};



#endif /* REAL_SLICER_OV_H_ */

