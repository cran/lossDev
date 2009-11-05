#ifndef REAL_SLICER_OV_H_
#define REAL_SLICER_OV_H_

#include <BASE/samplers/RealSlicer.h>

class StochasticNode;


/**
 * Slice sampler for real-valued distributions
 */
class RealSlicerOV : public base::RealSlicer
{
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
    RealSlicerOV();
    //virtual double value() const;
    //virtual void setValue(double value);
    //virtual void getLimits(double *lower, double *upper) const;
    virtual void update(RNG *rng);
    virtual void updateOverrelaxed(RNG *rng);
    virtual void localUpdateStep(RNG *rng);
    virtual std::string name() const;
    static bool canSample(StochasticNode const *node);
    //virtual bool isAdaptive();
    virtual bool adaptOff();
};



#endif /* REAL_SLICER_OV_H_ */

