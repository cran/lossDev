/*
##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, National Council On Compensation Insurance Inc.,                        ##
##                                                                                              ##
##    This file is part of lossDev.                                                             ##
##                                                                                              ##
##    lossDev is free software: you can redistribute it and/or modify                           ##
##    it under the terms of the GNU General Public License as published by                      ##
##    the Free Software Foundation, either version 3 of the License, or                         ##
##    (at your option) any later version.                                                       ##
##                                                                                              ##
##    This program is distributed in the hope that it will be useful,                           ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of                            ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             ##
##    GNU General Public License for more details.                                              ##
##                                                                                              ##
##    You should have received a copy of the GNU General Public License                         ##
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.                     ##
##                                                                                              ##
##################################################################################################
*/

#include "RJumpSplineFactory.h"
#include "RJumpSpline.h"
#include <JAGS/graph/Node.h>
#include <JAGS/graph/ConstantNode.h>
#include <JAGS/sampler/Linear.h>
#include <vector>
#include <set>

using std::set;
using std::vector;

RJumpSplineFactory::RJumpSplineFactory()
{
}

RJumpSplineFactory::~RJumpSplineFactory()
{
}

//needs work
bool RJumpSplineFactory::canSample(StochasticNode *node, Graph const &graph) 
const
{
	if( node->distribution()->name() != "dspline")
		return false;
	vector<Node const*>::const_iterator p;
	vector<Node const*> const &par = node->parents();
	

	
	//this isn't working
	for(p = par.begin(); p != par.end(); ++p)
	{
		Node const &n = **p;
		if(!static_cast<ConstantNode const *>(&n))
			return false;
	}
	
	
    vector<StochasticNode const*> stoch_nodes;
    vector<Node*> dtrm_nodes;
    Sampler::classifyChildren(vector<StochasticNode*>(1,node), 
		              graph, stoch_nodes, dtrm_nodes);
	

    // Check stochastic children
    for (unsigned int i = 0; i < stoch_nodes.size(); ++i) 
    {
    	StochasticNode const &sn = *stoch_nodes[i];
    	if(sn.distribution()->name() != "dnorm")
    		return false;
    	
    	if (isBounded(&sn))
    		return false; //Truncated distribution
    }

    // Check linearity of deterministic descendants
    if (!checkLinear(vector<StochasticNode*>(1, node), graph, false))
    	return false;


    
	return true;
}

Sampler * RJumpSplineFactory::makeSingletonSampler(StochasticNode *node,
				  Graph const &graph) const
{
	std::vector<StochasticNode*> vnode;
	vnode.push_back(node);	
	return  new RJumpSpline(vnode, graph);
}

