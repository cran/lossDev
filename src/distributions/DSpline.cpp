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

#include "DSpline.h"

#include <stdexcept>
#include <vector>

#include <util/nainf.h>
#include <JRmath.h>

#include <cmath>

using std::vector;

#include <iostream>
DSpline::DSpline():ArrayDist("dspline", 11)
{
	
}

DSpline::~DSpline()
{

}

double 
DSpline::logLikelihood(double const *x, unsigned int length,
                       vector<double const *> const &parameters,
                       vector<vector<unsigned int> > const &dims,
                       double const *lbound, double const *ubound) const
{
  //throw std::logic_error("ll doesn't make sense for dspline.  Neither dspline nor any node around dspline can use a default sampler.");
  //only cacl ll for 3rd and 4th error terms so we can slice sample their precisions (7th and 8th parameters) all other uses of ll will be "INVALID"

  const unsigned int &n = length;
  const unsigned int numberOfSplines =  dims[1][0];
  const unsigned int ncol = numberOfSplines == 1 ? 4 : 6;
  const unsigned int nrow = n / ncol;

  const double &tau1 = parameters[6][0];
  const double &tau2 = parameters[8][0];

  const double &rho = parameters[9][0];
  const double &etaRho = parameters[10][0];

  double ans = 0.0;

  if(numberOfSplines == 1)
  {
      //3rd column of error terms (prec for first value doesn't matter)
      //first value is head of AR1
      ans += dnorm4(x[1 + nrow * 2], 0.0, std::sqrt(1.0/tau1 / (1.0 - etaRho * etaRho)), 1 );
      for(unsigned int i = 2; i < nrow; ++i)
      	  ans += dnorm4(x[i + nrow * 2], etaRho * x[i-1 + nrow*2] ,std::sqrt(1.0/tau1), 1);
      
      
      //4th column of error terms (first value is blank)
      ans += dnorm(x[1 + nrow * 3], 0.0, std::sqrt(1.0/tau2 / (1.0 - rho * rho)), 1);
      for(unsigned int i = 2; i < nrow; ++i)
        ans += dnorm4(x[i + nrow * 3], rho * x[i-1 + nrow*3],std::sqrt(1.0/tau2), 1);
    }
  else if(numberOfSplines == 2)
    {
	//3rd column of error terms (prec for first value doesn't matter)
	//first value is head of AR1
	ans += dnorm4(x[1 + nrow * 4], 0.0, std::sqrt(1.0/tau1 / (1.0 - etaRho * etaRho)), 1 );
	for(unsigned int i = 2; i < nrow; ++i)
	    ans += dnorm4(x[i + nrow * 4], etaRho * x[i-1 + nrow*4] ,std::sqrt(1.0/tau1), 1);
	
	//4th column of error terms (first value is blank)
	ans += dnorm(x[1 + nrow * 5], 0.0, std::sqrt(1.0/tau2 / (1.0 - rho * rho)), 1);
	for(unsigned int i = 2; i < nrow; ++i)
	    ans += dnorm4(x[i + nrow * 5], rho * x[i-1 + nrow*5],std::sqrt(1.0/tau2), 1);
    }
  else
    {
      throw std::logic_error("Error in logLikeLihood for DSpline: numberOfSplines (:= dims[1][0]) must be equal to 1 or 2");
    }

  
  return ans;
    
    
}


vector<unsigned int> DSpline::dim (vector <vector<unsigned int> > const &args) const
{
  std::vector<unsigned int> dim;
  //number of rows equal to length of first parameter vector
  dim.push_back(args[0][0]);
  /*four columns in total
   *first is smoothed values
   *second is all set to equal the number of knots
   *third is a vector of normaly dist errors with zero mean and percision equal to the seventh parameter
   *fourth is a vector of normaly dist errors with zero mean and percision equal to the eighth parameter
   */
  dim.push_back(args[1][0] == 1 ? 4 : 6);
  return dim;
}

bool 
DSpline::checkParameterDim (vector<vector<unsigned int> > const &parameters) const
{

  //this should be checked by JAGS, see constructor
  //if (parameters.size() != 10)
  //return false;
	
  //std::cout << "starting dim test1" << std::endl;
  //first dim is number of evaluations
  if(parameters[0].size() != 1 || parameters[0][0] < 4)
    return false;

  //std::cout << "starting dim test2" << std::endl;
  //second is upper limit for number of knots should be one or two
  if(parameters[1].size() != 1 || (parameters[1][0] != 1 && parameters[1][0] != 2))
    return false;

  unsigned int numberOfSplines = parameters[1][0];

  //std::cout << "starting dim test3" << std::endl;
  //third is parm for tau
  if(parameters[2].size() != 1 || parameters[2][0] != 1)
    return false;

  //std::cout << "starting dim test4" << std::endl;
  //fourth is lower limit for knots values
  if(parameters[3].size() != 1 || parameters[3][0] != numberOfSplines)
    return false;

  //std::cout << "starting dim test5" << std::endl;
  //fifth is upper limit for knot values
  if(parameters[4].size() != 1 || parameters[4][0] != numberOfSplines)
    return false;

  //std::cout << "starting dim test6" << std::endl;
  //sixth is a and b for beta prior on knot positions
  if(numberOfSplines == 2)
    {
      if(parameters[5].size() != 2 || parameters[5][0] != 2 || parameters[5][1] != numberOfSplines)
        return false;
    }
  else if(numberOfSplines == 1)
    {
      if(parameters[5].size() != 1 || parameters[5][0] != 2)
        return false;
    }
  else
    {
      return false;
    }

    //std::cout << "starting dim test7" << std::endl;
  //seventh is the precision for the 3rd column which is a vector of normally dist error terms with zero mean first value gets its own nonstochastic precision
  if(parameters[6].size() != 1 || parameters[6][0] != 1)
    return false;

  //std::cout << "starting dim test8" << std::endl;
  //eighth is the precision for the first value in the 4th column which must be nonstochastic
  if(parameters[7].size() != 1 || parameters[7][0] != 1)
    return false;

  //std::cout << "starting dim test9" << std::endl;
  //nineth is the precision for the 4th column which is a vector of normally dist error terms with zero mean first value is blank
  if(parameters[8].size() != 1 || parameters[8][0] != 1)
    return false;  

  
  //tenth is the rho for the 4th column;
  if(parameters[9].size() != 1 || parameters[9][0] != 1)
    return false; 
 
  //eleventh is the rho for the 3rd column;
  if(parameters[10].size() != 1 || parameters[10][0] != 1)
    return false; 


  //std::cout << "dim test OK" << std::endl;
  return true;
}


bool 
DSpline::checkParameterValue(vector<double const *> const &parameters,
                             vector<vector<unsigned int> > const &dims) const
{
	

    //unsigned int numberOfSplines = parameters[1][0] == 1 ? 4 : 6;

  //std::cout << "starting test1" << std::endl;
  //x vector
  for(unsigned int i = 1; i < dims[0].size(); ++i)
    if(parameters[0][i-1] >= parameters[0][i])  //make sure x is sorted
      return false;


  //std::cout << "starting test2" << std::endl;
  //upper bound for number of knots
  for(unsigned int i = 0; i < dims[1][0]; ++i)
    {
      if(parameters[1][i] < 0  ||  //must not have negative number of knots 
         parameters[1][i] != static_cast<int>(parameters[1][i])) //upper bound must be an int
        return false;
    }

  //std::cout << "starting test3" << std::endl;
  //parms for tau
  if(parameters[2][0] <= 0) //make sure tau is greater than zero
    return false;
  
  //std::cout << "starting test4" << std::endl;
  //limit for knot values
  for(unsigned int i = 0; i < dims[3][0]; ++i)
    {
      if(parameters[4][i] <= parameters[3][i])
        return false;
    }
  
  
  
  //std::cout << "starting test5" << std::endl;
  // a and b for beta prior on knot positions

  unsigned int length = dims[5].size() == 2 ? dims[5][0] * dims[5][1] : dims[5][0];
  for(unsigned int i = 0; i < length; ++i)
    {
      if(parameters[5][i] <= 0 )
        return false;
    }

    

  //std::cout << "starting test6" << std::endl;
  //parms for tau
  if(parameters[6][0] <= 0) //make sure tau is greater than zero
    return false;

  //std::cout << "starting test7" << std::endl;
  //parms for tau
  if(parameters[7][0] <= 0) //make sure tau is greater than zero
    return false;

  //std::cout << "starting test8" << std::endl;
  //parms for tau
  if(parameters[8][0] <= 0) //make sure tau is greater than zero
    return false;

  if(parameters[9][0] <= -1 ||parameters[9][0] >= 1) //make sure rho produces a stationary result
      return false;

  if(parameters[10][0] <= -1 ||parameters[10][0] >= 1) //make sure rho produces a stationary result
      return false;

  //std::cout << "no problem" << std::endl;
  return true;
}

void 
DSpline::randomSample(double *x, unsigned int length,
                      std::vector<double const *> const &parameters,
                      std::vector<std::vector<unsigned int> > const  &dims,
                      double const *lbound, double const *ubound, RNG *rng) const
{
  //typicalValue(x, length, parameters, dims, lbound, ubound);
  throw std::logic_error("Doesn't make sense to random sample from a spline\nDid you remember to set initial values for dspline?");
}

void 
DSpline::typicalValue(double *x, unsigned int length,
                      std::vector<double const *> const &parameters,
                      std::vector<std::vector<unsigned int> > const &dims,
                      double const *lbound, double const *ubound) const
{

  //vector<unsigned int> const &d = dim(dims);
	
  /*
    for(unsigned int i = 0; i < d[0]; ++i)
    x[i] = 0;
	
    for(unsigned int i = d[0]; i < 2*d[0]; ++i)
    x[i] = 0;
  */
	
  for(unsigned int i = 0; i < length; ++i)
    x[i] = 0;
	


}


void DSpline::support(double *lower, double *upper, unsigned int length,
                      std::vector<double const *> const &parameters,
                      std::vector<std::vector<unsigned int> > const &dims) const
{
  for(unsigned int i = 0; i < length; ++i)
    {
      lower[i] = JAGS_NEGINF;
      upper[i] = JAGS_POSINF;
    }
	
}

bool DSpline::isSupportFixed(std::vector<bool> const &fixmask) const
{
  return true;
}





