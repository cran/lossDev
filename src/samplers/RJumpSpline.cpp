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

#include <graph/Node.h>
#include <graph/StochasticNode.h>
#include <graph/DeterministicNode.h>
#include <graph/StochasticNode.h>
#include <graph/LogicalNode.h>
#include <distribution/Distribution.h>
#include <function/Function.h>
#include <graph/NodeError.h>
#include <RNG.h>
#include <sampler/GraphView.h>
#include <module/ModuleError.h>

#include "MNorm.h"

#include <cmath>

#include "RJumpSpline.h"

#include <JRmath.h>

#include <vector>
#include <set>
#include <iostream>
#include <cstring>

#include <R_ext/Lapack.h>

using std::vector;
using std::set;
using std::memcpy;

//calcBeta for 
//this is taken from JAGS/bus/MNormal
void RJumpSpline::calBeta(double *betas, bool const &current, unsigned int const &chain)
{

  
  double *coeff;	
	
  if(current)
    {
      coeff = _currentBeta[chain];
		
    } else
    {
      coeff = _proposedBeta;
    }
  //for efficiency assume that the value is set to whatever value current indicates
  //This is safe since calPost is the only function to call calBeta and it sets the value before it calls calBeta
  //We don't know for sure the value of the node is set to the value we are talking about
  //setSplineValue(*K,xnew, T, chain);
	
  //double const *xold = coeff;
  unsigned int nrow = betaLength(current, chain);
  //unsigned int nrow = snode->length();

  /*double *xnew = new double[nrow];
  for (unsigned int i = 0; i < nrow; ++i) 
    {
      xnew[i] = xold[i];
    }
  */

  vector<StochasticNode const*> const &stoch_children = 
      _mgv->stochasticChildren();
  
  unsigned long nchildren = stoch_children.size();
  double *beta_j = betas;

  for (unsigned int j = 0; j < nchildren; ++j) 
    {
      StochasticNode const *snode = stoch_children[j];
      double const mu = snode->parents()[0]->value(chain)[0];
      
      for (unsigned int i = 0; i < nrow; ++i) 
      {
	  beta_j[i] = -mu;
      }
      
      beta_j += nrow;
    }

    
  for (unsigned int i = 0; i < nrow; ++i) 
    {
      coeff[i] += 1;
      //sampler->setValue(xnew, nrow, chain);
      setSplineValue(current, chain);
      beta_j = betas;
      for (unsigned int j = 0; j < nchildren; ++j) 
        {
          StochasticNode const *snode = stoch_children[j];
          double const mu = snode->parents()[0]->value(chain)[0];
	  beta_j[i] += mu;
          beta_j += nrow;
        }
      coeff[i] -= 1;
    }
  //sampler->setValue(xnew, nrow, chain);
    
  setSplineValue(current, chain);
  //delete [] xnew;
}



//this is taken from JAGS/bugs/MNormal

void RJumpSpline::calPost(bool const &current, unsigned int chain)
{
  double const *coeff;
  double * b;
  double * A;
	
  if(current)
    {
      coeff = _currentBeta[chain];
      b = _bPostCurrent;
      A = _APostCurrent;
		
    } else
    {
      coeff = _proposedBeta;
      b = _bPostProposed;
      A = _APostProposed;
    }	

  //We don't know for sure the value of the node is set to the value we are talking about
  setSplineValue(current, chain);
  
  vector<StochasticNode const*> const &stoch_children = _mgv->stochasticChildren();
  unsigned int nchildren = stoch_children.size();
  
  unsigned int length_betas = nchildren * betaLength(current, chain);
  
  StochasticNode const *snode = _snode;
  double const *xold = coeff;

  int nrow =   betaLength(current, chain);
  unsigned int nrowUS = betaLength(current, chain);
  /* 
     The log of the full conditional density takes the form
     -1/2(t(x) %*% A %*% x - 2 * b %*% x)
     
     For computational convenience, we reset the origin to xold,
     the current value of the node.
  */
  int N = nrow * nrow;


  for(unsigned int i = 0; i < nrowUS; ++i)
  {
      if(i == 0)
	  b[i] =   _tauOfFirstIn3rdColumn * (0.0 - xold[i]);
      else if(i < 2) // we are in the first _n - 1 values (which are the error in the 3rd column)
	  b[i] =  _tauOf3rdColumn * (1 - _etaRho * _etaRho) * (0.0 - xold[i]);
      else if(i < 1 + (_TriDim - 1))
	  b[i] = _tauOf3rdColumn * (0.0 - xold[i]);
      else if(i < 1 + 1 + (_TriDim-1)) // we are in the second _n - 1 values (which are the errors in the 4th column)
	  b[i] = _tauOf4thColumn * (1 - _rho * _rho) * (0.0 - xold[i]);
      else if(i < 1 + 2 * (_TriDim-1))
        b[i] = _tauOf4thColumn * (0.0 - xold[i]);
      else  //we are talking about an intercept slope or spline coeff
        b[i] = _tau * (0.0 - xold[i]);
    }

  for(unsigned int i = 0; i < nrowUS; ++i)
    for(unsigned int j = 0; j < nrowUS; ++j)
      {
        if(i == j) // we are on the diagional
          {
            if(i == 0)
              A[i + j * nrow] =   _tauOfFirstIn3rdColumn;
	    else if (i < 2)
		A[i + j * nrow] =  _tauOf3rdColumn * (1 - _etaRho * _etaRho);
            else if(i < 1 + (_TriDim - 1)) // we are in the first _n - 1 values (which are the error in the 3rd column)
              A[i + j * nrow] =  _tauOf3rdColumn;
	    else if(i < 1 + 1 + (_TriDim-1)) // we are in the second _n - 1 values (which are the errors in the 4th column)
		A[i + j * nrow]	 = _tauOf4thColumn * (1 - _rho * _rho) ;
            else if(i < 1 + 2 * (_TriDim-1)) // we are in the second _n - 1 values (which are the errors in the 4th column)
              A[i + j * nrow] = _tauOf4thColumn;
            else  //we are talking about an intercept slope or spline coeff
              A[i + j * nrow] = _tau ;
            
          }
        else  // on the off diagional
          A[i + j * nrow] = 0;
      }
  
  
  /* FORTRAN routines are all call-by-reference, so we need to create
   * these constants */
   int i1 = 1;

  _calPost_betas->makeSufficientWithoutCopy(length_betas);
  double* betas=_calPost_betas->value(); //betas = new double[length_betas];
  calBeta(betas, current, chain);
  
  //Calculate largest possible size of working matrix C
  //Chris Laws Changed this from JAGS because it looked like a bug
  //My answer will be >= original so should be safe
  //Need to double check
  
  /* Now add the contribution of each term to A, b 
     
  b += N_j * beta_j %*% tau_j (Y_j - mu_j)
  A += N_j * beta_j %*% tau_j %*% t(beta_j)
  
  where 
  - N_j is the frequency weight of child j
  - beta_j is a matrix of linear coefficients
  - tau_j is the variance-covariance matrix of child j
  - mu_j is the mean of child j
  - Y_j is the value of child j
  
  We make use of BLAS routines for efficiency.
  
  */
  double const *beta_j = betas;
  for (unsigned int j = 0; j < nchildren; ++j) 
    {
      
      StochasticNode const *snode = stoch_children[j];
      double const *Y = snode->value(chain);
      double const *mu = snode->parents()[0]->value(chain);
      double const *tau = snode->parents()[1]->value(chain);
   
      
      
      double alpha = tau[0];
      F77_NAME(dsyr)("L", &nrow, &alpha, beta_j, &i1, A, &nrow);
      alpha *= (Y[0] - mu[0]);
      F77_NAME(daxpy)(&nrow, &alpha, beta_j, &i1, b, &i1);
        

      beta_j += nrow;
      
    }
 
  
  
  /* 
     Solve the equation A %*% x = b to get the posterior mean.
     We have to take a copy of A as it is overwritten during
     the call to DPOSV. The result is stored in b
  */

  _calPost_Acopy->makeSufficientWithoutCopy(N);
  double* Acopy=_calPost_Acopy->value();
  memcpy(Acopy, A, sizeof(double) * N);

  int one = 1;
  int info;
  F77_NAME(dposv) ("L", &nrow, &one, Acopy, &nrow, b, &nrow, &info);
  if (info != 0) 
    {
      throwNodeError(snode,
                      "unable to solve linear equations in Conjugate mnorm sampler");
    }
  
  //Shift origin back to original scale
  for (int i = 0; i < nrow; ++i) 
    b[i] += xold[i];
  	
}




RJumpSpline::RJumpSpline(GraphView *gv):
    Sampler(gv),_mgv(gv)
{
    _snode = gv->nodes()[0];
    _nchain = _snode->nchain();
    _numberOfSplines = _snode->parents()[1]->dim()[0];
    
    _knots= new Knots*[_numberOfSplines];
    
    double PriorForT[2];
    //unsigned int KLimits[2];
    double TLimits[2];
    
    vector<unsigned int> maxK(_numberOfSplines);
    for(unsigned int i = 0; i < _numberOfSplines; ++i)
    {
	
	PriorForT[0] = _snode->parameters(0)[5][0 + 2 * i];
	PriorForT[1] = _snode->parameters(0)[5][1 + 2 * i];
	
	//KLimits[0] = 0;
	//KLimits[1] = static_cast<unsigned int>(_snode->parameters(0)[1][i]);

	TLimits[0] = _snode->parameters(0)[3][i];
	TLimits[1] = _snode->parameters(0)[4][i];
	
	
	_knots[i] = new Knots(_nchain, //number of chains
			      PriorForT,//prior for knot positions
			      gv->nodes()[i + 1],//number of knots
			      TLimits//Limits for knot positions
	    );
	
	maxK[i] = _knots[i]->maxK();
    }
    
    
    //std::cout <<"mid Const1" << std::endl;
    

    //values at which spline is to be evaluated
    //For now assume this is constant
    
    _n = new unsigned int[_numberOfSplines];
    
    for(unsigned int i = 0; i < _numberOfSplines; ++i)
	_n[i] = static_cast<unsigned int>(_snode->parents()[0]->dim()[0]);
    
    
    _TriDim = _n[0];
    
    _x = new double*[_numberOfSplines];
    
    for(unsigned int i = 0; i < _numberOfSplines; ++i)
    {
	_x[i] = new double[_n[i]];
	for(unsigned int j = 1; j <= _n[i]; ++j)
	    _x[i][j-1] = j;
    }
    
    //std::cout <<"mid Const2" << std::endl;		
    //Coefficents
    //first value is an error term with zero mean and precision defined by the 8th parameter with must be NON-STOCHASTIC (_tauOfFirstIn3rdColumn)
    //next _TriDim-1 values are for the 3rd column of error terms after the first element
    //next _TriDim-1 values are for the 4th column of error terms (first value of this column is blank)
    //next value is intercept
    //next value is slope
    //next value is intercept
    //next value is slope
    // length of Coefficents is thus = 1 + 2 *(_TriDim-1) + 2 * _numberOfSplines + number of knots {+ number of knots}
    _currentBeta = new double*[_nchain];
    unsigned int maxBetaLength = 1 + 2*(_TriDim-1) + 2 * _numberOfSplines;
    for(unsigned int i = 0; i < _numberOfSplines; ++i)
	maxBetaLength += maxK[i];
    
  
    for(unsigned int i = 0; i < _nchain; ++i)
    {
	
	_currentBeta[i] = new double[maxBetaLength];
	for(unsigned int j = 0; j < maxBetaLength; ++j)
	    _currentBeta[i][j] = 0;
    }
    
    _proposedBeta = new double[maxBetaLength];
    
    //Prior for Betas and intercept
    _tau = _snode->parameters(0)[2][0];
    
    //Prior precision for 3rd and 4th column error terms;
    //initially set these to be the value in the first chain, update will have to reset for each chain;
    _tauOf3rdColumn = _snode->parameters(0)[6][0];
    _tauOf4thColumn = _snode->parameters(0)[8][0];
    
    //Prior precision for first value in the 3rd column NON-STOCHASTIC
    _tauOfFirstIn3rdColumn= _snode->parameters(0)[7][0];

    //autoregressive coefficient
    _rho =  _snode->parameters(0)[9][0];
    _etaRho =  _snode->parameters(0)[10][0];
    
    //posterior parameters for Beta under current and posterior cases
    _bPostCurrent = new double[maxBetaLength];
    _APostCurrent = new double[maxBetaLength * maxBetaLength];
  
    _bPostProposed = new double[maxBetaLength];
    _APostProposed = new double[maxBetaLength * maxBetaLength];
    
    _updatingKnotsI = 0;
    for(unsigned int i = 0; i < _nchain; ++i)
	setSplineValue(true, i);
    
    
    _calPost_betas = new ExpandableArray;
    _calPost_Acopy = new ExpandableArray;
    
}

RJumpSpline::~RJumpSpline()
{
    
   for(unsigned int i = 0; i < _numberOfSplines; ++i)
    delete _knots[i];
   delete[] _knots;
  
  delete[] _n;

  for(unsigned int i = 0; i < _numberOfSplines; ++i)
    delete[] _x[i];
  
  delete[] _x;
  
	
  //11/18/08 Chris Laws added following two lines to avoid memory leaks
  //shouldn't make much of a difference though because these are only called once jags is done
  for(unsigned int i = 0; i < _nchain; ++i)
    delete[] _currentBeta[i];
  delete[] _currentBeta;
  delete[] _proposedBeta;
	


  //11/18/08 Chris Laws added following two lines to avoid memory leaks
  //shouldn't make much of a difference though because these are only called once jags is done	
  delete [] _bPostCurrent;
  delete [] _APostCurrent;
	
  delete [] _bPostProposed;
  delete [] _APostProposed;
	
  delete _calPost_betas;
  delete _calPost_Acopy;
}

double RJumpSpline::llZ(bool const &current, unsigned int chain) const
{

  //std::cout << "start of llZ" << std::endl;
  double const *coeff;
  double const *b;
  double const *A;
  
  if(current)
    {
      coeff = _currentBeta[chain];
      b = _bPostCurrent;
      A = _APostCurrent;
      
    } else
    {
      coeff = _proposedBeta;
      b = _bPostProposed;
      A = _APostProposed;
    }
  //We don't know for sure the value of the node is set to the value we are talking about
  setSplineValue(current,chain);

  double ans = 0.0;
   unsigned int length = betaLength(current, chain);
  for(unsigned int i = 0; i < length; ++i)
    {
      if(i == 0)
        ans  +=  dnorm4(coeff[i], 0.0, std::sqrt(1.0/ _tauOfFirstIn3rdColumn), 1);
      
      else if(i < 2)
	  ans  +=  dnorm4(coeff[i], 0.0, std::sqrt(1.0/ _tauOf3rdColumn / (1 - _etaRho * _etaRho)), 1);
      else if(i < 1 + (_TriDim - 1)) // we are in the first _n - 1 values (which are the error in the 3rd column)
	  ans  +=  dnorm4(coeff[i], 0.0, std::sqrt(1.0/ _tauOf3rdColumn), 1);
      
      else if(i < 1 + 1 + (_TriDim-1)) // we are in the second _n - 1 values (which are the errors in the 4th column)
	  ans += dnorm4(coeff[i], 0.0, std::sqrt(1.0 / _tauOf4thColumn / (1 - _rho * _rho)), 1) ;

      else if(i < 1 + 2 * (_TriDim-1)) // we are in the second _n - 1 values (which are the errors in the 4th column)
        ans  +=  dnorm4(coeff[i], 0.0, std::sqrt(1.0/ _tauOf4thColumn), 1);
      
      else  //we are talking about an intercept slope or spline coeff
        ans  +=  dnorm4(coeff[i], 0.0, std::sqrt(1.0/ _tau), 1);
    }
  
  vector<StochasticNode const *> const & children = _mgv->stochasticChildren();
  
  for(unsigned int i = 0; i < children.size(); ++i)
    {
      StochasticNode const &sn = *children[i];
      ans += sn.logDensity(chain, PDF_FULL);
    }
  
  
  vector<double const*> par;
  par.push_back(b);
  par.push_back(A);
  
  vector<vector<unsigned int> > d;


  vector<unsigned int> m;
  m.push_back(length);
  d.push_back(m);
  
  vector<unsigned int> per;
  per.push_back(length);
  per.push_back(length);
  d.push_back(per);
  
  //assume the coeff are areally set
  ans -= MNorm_logLikelihood(coeff, length,
			     par,
			     d,
			     NULL,NULL);
  
  
  //setSplineValue(*K, coeff, T, chain);
  
  //std::cout << "end of llZ" << std::endl;
  return ans;
}



void RJumpSpline::accept(unsigned int const &chain, unsigned int const &spline, Knots::TypeOfUpdate type, RNG * const & rng)
{

  double p = llZ(false, chain) - llZ(true, chain);
  p += _knots[spline]->acceptProbBalance(chain, type);
    
  if(rng->uniform() < std::exp(p))
    {
      unsigned int length = betaLength(false, chain);
      for(unsigned int i = 0; i < length; ++i)
        _currentBeta[chain][i] = _proposedBeta[i];

      _knots[_updatingKnotsI] -> acceptProposedValues(chain);
    }
  
  setSplineValue(true, chain);
	
}

void RJumpSpline::update(std::vector<RNG *> const &rng)
{
  
  //std::cout << "start of update" << std::endl;
  for(unsigned int c = 0; c < _nchain; ++c)
    {
      
      //Prior precision for 3rd and 4th column error terms;
      _tauOf3rdColumn = _snode->parameters(c)[6][0];
      _tauOf4thColumn = _snode->parameters(c)[8][0];

      //autoregressive coefficient
      _rho =  _snode->parameters(c)[9][0];
      _etaRho = _snode->parameters(c)[10][0];
        
      RNG* const r = rng[c];

      for(unsigned int i = 0; i < _numberOfSplines; ++i)
        {
          _updatingKnotsI = i;
          Knots::TypeOfUpdate type = _knots[i]->createProposal(r, c);

    

          if(type == Knots::Nothing)
            {
              calPost(true, c);
              MNorm_randomsample(_currentBeta[c], _bPostCurrent, _APostCurrent, true, betaLength(true, c), r);
              setSplineValue(true, c);
              continue;
            }

          calPost(true, c);
          MNorm_randomsample(_currentBeta[c], _bPostCurrent, _APostCurrent, true, betaLength(true, c), r);
          calPost(false, c);
          MNorm_randomsample(_proposedBeta, _bPostProposed, _APostProposed, true, betaLength(false, c), r);
          accept(c, i, type, r);
        }
    }
  //std::cout << "end of update\n" << std::endl;
	
}

bool RJumpSpline::checkAdaptation() const
{
  return true;
}



void RJumpSpline::adaptOff()
{

}

bool RJumpSpline::isAdaptive() const
{
  return false;
}



std::string RJumpSpline::name() const
{
  return "RJumpSpline";
}



void RJumpSpline::setSplineValue(bool const &current, unsigned int const &chain) const
{


  double const *coeff;

  if(current)
    coeff = _currentBeta[chain];
  else
    coeff = _proposedBeta;
    
  unsigned int numberOfKnots;
  double tmpKnotPos;
  
  const unsigned int &nrow = _TriDim;
  const unsigned int ncol = _numberOfSplines == 1 ? 4 : 6;

  double * value = new double[nrow * ncol + _numberOfSplines];

  //first spline
  if(current || _updatingKnotsI != 0)
    {

      numberOfKnots = _knots[0]->currentLength(chain);

    }
  else
    numberOfKnots = _knots[0]->proposedLength(chain);

  
  for(unsigned int i = 0; i < _n[0]; ++i)
    {

      value[i + nrow * 0] = coeff[0 + 1 + 2*(_TriDim-1)] + coeff[1 + 1 + 2*(_TriDim-1)] * (_x[0][i] - _knots[0]->minT());
      
      for(unsigned int j = 0; j < numberOfKnots; ++j)
        {

          if(current || _updatingKnotsI != 0)
            tmpKnotPos = _knots[0]->currentValue(chain, j);
          else
            tmpKnotPos = _knots[0]->proposedValue(chain, j);
          
          if(_x[0][i] - tmpKnotPos <= 0)
            break;
          //don't forget to scale up the knot positions
          value[i + nrow * 0] += coeff[j + 1 + 2*(_TriDim - 1) + 2 * _numberOfSplines] * (_x[0][i] - tmpKnotPos);  
          
        }
    }
  unsigned int offset = numberOfKnots;
  
  //second column number of knots
  for(unsigned int i = 0; i < nrow; ++i)
    value[i + nrow * 1] = numberOfKnots;
  value[nrow * ncol] = numberOfKnots;

  if(_numberOfSplines == 2)
  {


      //second spline
      if(current || _updatingKnotsI != 1)
        numberOfKnots = _knots[1]->currentLength(chain);
      else
        numberOfKnots = _knots[1]->proposedLength(chain);
  
      for(unsigned int i = 0; i < _n[1]; ++i)
        {
          value[i + nrow * 2] = coeff[2 + 1 + 2*(_TriDim-1)] + coeff[3 + 1 + 2*(_TriDim-1)] * (_x[1][i] - _knots[1]->minT());
      
          for(unsigned int j = 0; j < numberOfKnots; ++j)
            {
              if(current || _updatingKnotsI != 1)
                tmpKnotPos = _knots[1]->currentValue(chain, j);
              else
                tmpKnotPos = _knots[1]->proposedValue(chain, j);
          
              if(_x[1][i] - tmpKnotPos <= 0)
                break;
              //don't forget to scale up the knot positions
              value[i + nrow * 2] += coeff[j + 1 + 2*(_TriDim - 1) + 2*2 + offset] * (_x[1][i] - tmpKnotPos);  
          
            }
        }
    


      //forth column number of knots
      for(unsigned int i = 0; i < nrow; ++i)
        value[i + nrow * 3] = numberOfKnots;
      value[nrow * ncol + 1] = numberOfKnots;

      //fifth column error terms first value blank
      value[0 + nrow * 4] = coeff[0]; 
      value[1 + nrow * 4] = coeff[1];
      for(unsigned int i = 1; i < _TriDim-1; ++i)
        value[i + nrow * 4 + 1] = coeff[i + 1] + _etaRho * value[i + nrow * 4 + 1 - 1];

      //sixth column error terms first value blank
      value[0 + nrow * 5] = 0; 
      value[0 + nrow * 5 + 1] = coeff[0 + 1 + (_TriDim-1)];
      for(unsigned int i = 1; i < _TriDim-1; ++i)
	  value[i + nrow * 5 + 1] = coeff[i + 1 + (_TriDim-1)] + _rho * value[(i-1) + nrow * 5 + 1];
    }
  else
  {
      
      //third column error terms first value blank
      value[0 + nrow * 2] = coeff[0]; 
      value[1 + nrow * 2] = coeff[1];
      for(unsigned int i = 1; i < _TriDim-1; ++i)
        value[i + nrow * 2 + 1] = coeff[i + 1] + _etaRho * value[i + nrow * 2 + 1 - 1];

      //fourth column error terms first value blank
      value[0 + nrow * 3] = 0; 
      value[0 + nrow * 3 + 1] = coeff[0 + 1 + (_TriDim-1)];
      for(unsigned int i = 1; i < _TriDim-1; ++i)
	  value[i + nrow * 3 + 1] = coeff[i + 1 + (_TriDim-1)] + _rho *  value[(i-1) + nrow * 3 + 1];
  }
    
  
  const_cast<RJumpSpline* const>(this)->_mgv->setValue(value, nrow * ncol + _numberOfSplines, chain);

   //std::cout << "end setSplineValue" << std::endl;

  delete [] value;
}

unsigned int RJumpSpline::betaLength(bool const &current, unsigned int const& chain) const
{
  unsigned int ans = 1 + 2 * (_TriDim - 1) + 2 * _numberOfSplines;

  
  if(current)
    {
      for(unsigned int i = 0; i < _numberOfSplines; ++i)
        ans  +=  _knots[i]->currentLength(chain);
    }
  else if (_updatingKnotsI == 0)
    {
      if ( _numberOfSplines == 2)
        ans += _knots[0]->proposedLength(chain) + _knots[1]->currentLength(chain);
      else
        ans += _knots[0]->proposedLength(chain);
    }  
  //this can only happen if we have two splines
  else if (_updatingKnotsI == 1)
    {
      ans += _knots[0]->currentLength(chain) + _knots[1]->proposedLength(chain);
    }

  return ans;
  
}
