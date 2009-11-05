//#include <config.h>

#include "RealSlicerOV.h"
//#include "DiscreteSlicer.h"
#include "SliceFactoryOV.h"

#include <sampler/DensitySampler.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>

#include <vector>

using std::vector;



bool 
SliceFactoryOV::canSample(StochasticNode * node, Graph const &graph) const
{
    return RealSlicerOV::canSample(node);
}

Sampler *SliceFactoryOV::makeSingletonSampler(StochasticNode *snode,
					      Graph const &graph) const
{
    unsigned int nchain = snode->nchain();
    vector<DensityMethod*> methods(nchain, 0);
    
    for (unsigned int ch = 0; ch < nchain; ++ch) {

	methods[ch] = new RealSlicerOV();

    }
    
    vector<StochasticNode*> nodes(1, snode);
    return new DensitySampler(nodes, graph, methods);
}


