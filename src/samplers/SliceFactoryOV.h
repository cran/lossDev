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
    bool canSample(StochasticNode *snode, Graph const &graph) const;
    Sampler *makeSingletonSampler(StochasticNode *snode, Graph const &graph)
	const;
};


#endif /* SLICE_FACTORY_OV_H_ */
