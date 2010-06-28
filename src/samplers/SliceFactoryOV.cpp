//#include <config.h>

#include "RealSlicerOV.h"
//#include "DiscreteSlicer.h"
#include "SliceFactoryOV.h"

#include <sampler/ParallelSampler.h>
#include <sampler/SampleMethod.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>

#include <vector>

using std::vector;



bool 
SliceFactoryOV::canSample(StochasticNode * node, Graph const &graph) const
{
    return RealSlicerOV::canSample(node);
}

Sampler *SliceFactoryOV::makeSampler(StochasticNode *snode,
				     Graph const &graph) const
{
    unsigned int nchain = snode->nchain();
    vector<SampleMethod*> methods(nchain, 0);

    vector<StochasticNode*> sv;
    sv.push_back(snode);
    GraphView * gv = new GraphView(sv, graph);
    
    for (unsigned int ch = 0; ch < nchain; ++ch) {

	methods[ch] = new RealSlicerOV(gv, ch);

    }
    
    vector<StochasticNode*> nodes(1, snode);
    return new ParallelSampler(gv,  methods);
}

std::string SliceFactoryOV::name() const
{
    return "SliceFactoryOV";
}



