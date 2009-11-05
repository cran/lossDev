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


#include "RJump.h"
#include <vector>

#include "distributions/DSpline.h"
//#include "distributions/DBlockTChisqr.h"
#include "distributions/DChisqrOV.h"
#include "distributions/DTOV.h"
#include "distributions/DNormOV.h"
#include "distributions/DUnifOV.h"


#include "samplers/RJumpSplineFactory.h"
//#include "samplers/MultiVarSliceFactory.h"
#include "samplers/SliceFactoryOV.h"



using std::vector;


RJump::RJump()
{


	insert(new DSpline);
	//insert(new DBlockTChisqr);
	insert(new DChisqrOV);
	insert(new DTOV);
	insert(new DNormOV);
	insert(new DUnifOV);


	insert(new RJumpSplineFactory);
	//insert(new MultiVarSliceFactory);
	insert(new SliceFactoryOV);
	
	
}

RJump::~RJump()
{
    vector<Distribution*> const &dvec = distributions();
    for (unsigned int i = 0; i < dvec.size(); ++i) 
    {
    	delete dvec[i];
    }
    



    vector<SamplerFactory*> const &svec = samplerFactories();
    for (unsigned int i = 0; i < svec.size(); ++i) {
	delete svec[i];
    }
}

RJump _RJump;


