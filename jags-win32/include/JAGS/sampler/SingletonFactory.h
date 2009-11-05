#ifndef SINGLETON_FACTORY_H_
#define SINGLETON_FACTORY_H_

#include <sampler/SamplerFactory.h>
/**
 * @short Factory object for a Sampler that samples a single node
 *
 * Many Sampler objects update a singe StochasticNode (and its
 * deterministic descendants).  This is a convenience class designed
 * to make it easier to create a factory object for such samplers.
 */
class SingletonFactory : public SamplerFactory
{
public:
    /**
     * Determines whether the factory can produce a Sampler for the
     * given node, within the given graph. This function is called
     * by SingletonFactory#makeSampler
     */
    virtual bool canSample(StochasticNode *node, Graph const &graph) 
	const = 0;
    /**
     * Returns a dynamically allocated Sampler for a given node. This
     * function is called by SingletonFactory#makeSampler.
     */
    virtual Sampler *makeSingletonSampler(StochasticNode *node,
					  Graph const &graph) const = 0;
    /**
     * This traverses the graph, creating a Sampler, when possible,
     * for each individual StochasticNode.
     */
    void makeSampler(std::set<StochasticNode*> &nodes, 
		     Graph const &graph,
		     std::vector<Sampler*> &samplers) const;
};

#endif /* SINGLETON_FACTORY_H */
