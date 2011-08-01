/*
##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, 2010, 2011 National Council On Compensation Insurance Inc.,             ##
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

#ifndef RJUMPSPLINEFACTORY_H_
#define RJUMPSPLINEFACTORY_H_

#include <vector>
#include <set>

#include <graph/Node.h>
#include <distribution/Distribution.h>
#include <graph/StochasticNode.h>
#include <sampler/SamplerFactory.h>


class RJumpSplineFactory : public SamplerFactory
{
public:
    RJumpSplineFactory();
    virtual ~RJumpSplineFactory();
    /**
     * Determines whether the factory can produce a Sampler for the
     * given node, within the given graph. This function is called
     * by SingletonFactory#makeSampler
     */
    //needs work
    virtual bool canSample(StochasticNode *node, Graph const &graph) 
	const;
    /**
     * Returns a dynamically allocated Sampler for a given node. This
     * function is called by SingletonFactory#makeSampler.
     */
    virtual std::vector<Sampler*> makeSamplers(std::set<StochasticNode*> const &nodes, 
					       Graph const &graph) const;

    virtual std::string name() const;
};

#endif /*RJUMPSPLINEFACTORY_H_*/
