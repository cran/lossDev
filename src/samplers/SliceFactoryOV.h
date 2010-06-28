#ifndef SLICE_FACTORY_OV_H_
#define SLICE_FACTORY_OV_H_

#include <sampler/SingletonFactory.h>
class StochasticNode;
class Graph;


/**
 * @short Factory object for slice samplers
 */
class SliceFactoryOV : public SingletonFactory
{
public:
    virtual bool canSample(StochasticNode *snode, Graph const &graph) const;
    virtual Sampler *makeSampler(StochasticNode *snode, Graph const &graph)
	const;
    virtual std::string name() const;
};


#endif /* SLICE_FACTORY_OV_H_ */
