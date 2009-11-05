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

#include <stdexcept>
#include <cmath>

#include <iostream>
#include "Knots.h"
#include <JAGS/JRmath.h>
#include <vector>
#include <set>

using std::vector;
using std::set;

static unsigned int factorial(unsigned int const &n)
{
  if(n <= 1)
    return 1;
	
  unsigned int ans = 1;
	
  for(unsigned int i = 2; i <= n; ++i)
    ans *= i;
	
  return ans;
		
		
}

unsigned int Knots::calcMaxDelta(unsigned int const &K, TypeOfUpdate type) const
{
	
  //enum TypeOfUpdate {Nothing, Birth, Death, Move};
  if(type == Nothing)
    return 0;
	
  if(type == Birth)
    return std::min(_maxK - K, static_cast<unsigned int>(7));
	
  if(type == Death)
    return std::min(K - _minK, static_cast<unsigned int>(7));
	
  if(type == Move)
    return K;
	
  throw std::logic_error("Tried to call calcMaxDelta With UnKnown TypeOfUpdate");
}


unsigned int Knots::rPoisT(double const &lambda, unsigned int const &T, RNG * const & rng) const
{
  if(T == 1)
    return 1;
	
  double p[T];
	
	
  p[0] = std::exp(- lambda);
	
  for(unsigned int i = 1; i < T; ++i)
    p[i] = p[i - 1] + std::exp(-lambda) * std::pow(lambda, static_cast<double>(i)) / factorial(i);
	
  for(unsigned int i = 0; i < T; ++i)
    p[i] /= p[T-1];
	
  double rn = rng->uniform();
	
  for(unsigned int i = 0; i < T; ++i)
    {
      if(rn < p[i])
        return i + 1;
    }
	
  throw std::runtime_error("not returning a valule in rPoisT");
		
}

double Knots::dPoisT(unsigned int const &x, double const &lambda, unsigned int const &T) const
{
	
  if(T == 1)
    {
      if(x == 1)
        return 0.0;
      else
        throw std::runtime_error("bad value in dPoisT");
    }
	

  double sum = std::exp(- lambda);
	
  for(unsigned int i = 1; i < T; ++i)
    sum += std::exp(- lambda) * std::pow(lambda, static_cast<double>(i)) / factorial(i);
	
  return std::log(std::exp(- lambda) * std::pow(lambda, static_cast<double>(T-1)) / factorial(T - 1) / sum);
	
	
}

double Knots::rBeta(double const &mu, double const &sig, RNG * const & rng) const
{
      
    double a = mu * pow(sig, -2);
    double b = (1.0 - mu) * pow(sig, -2);
    
    return rbeta(a, b, rng);
}


double Knots::dBeta(double const &x, double const &mu, double const &sig) const
{
      
    double a = mu * pow(sig, -2);
    double b = (1.0 - mu) * pow(sig, -2);
    
    return dbeta(x, a, b, 1);
}

void Knots::addToProposedT(double const &T)
{
  if(_proposedK == 0)
    {
      _proposedT[0] = T;
      ++_proposedK;
      return;
    }

  for(unsigned int i = _proposedK; i >=0; --i)
    {
		
      if(i == 0 || _proposedT[i-1] < T)
        {
          _proposedT[i] = T;
          ++_proposedK;
          break;
        }
      _proposedT[i] = _proposedT[i-1];
		
    }
}

double Knots::removeFromProposedT(unsigned int const &T)
{
  double ans = _proposedT[T];
  for(unsigned int i = T; i < _proposedK - 1; ++i)
    {
      _proposedT[i] = _proposedT[i+1];
    }
	
  --_proposedK;
	
  return ans;
}



Knots::Knots(unsigned int nchain, double const * PriorForT, unsigned int const * KLimits, double const * TLimits):
  _nchain(nchain),
  _aPriorForT(PriorForT[0]),  //a and b for beta prior of Knots
  _bPriorForT(PriorForT[1]),
  _minK(KLimits[0]), //min and max number of knots
  _maxK(KLimits[1]),
  _minT(TLimits[0]),    //min and max interval on which the Knots can fall //like in WinBUGS minT is also the position of Knot0 (which is not deemed to be part of _...T)	
  _maxT(TLimits[1])
  
{

  //std::cout << "beta priors are: " << _aPriorForT << " and " << _bPriorForT << std::endl;
  //number of knots
  _currentK = new unsigned int[_nchain];
  _proposedK = 0;
	
  for(unsigned int i = 0; i < _nchain; ++i)
    _currentK[i] = 1;
	
	
  //Position of Knots  these are stored on (0,1) and adjusted by _minT _maxT
  _currentT = new double*[_nchain];
  for(unsigned int i = 0; i < _nchain; ++i)
    {
      _currentT[i] = new double[_maxK];
      _currentT[i][0] = 0.5;
    }
	
  _proposedT = new double[_maxK];
	
  //Values needed to go from _currentT to _proposedT and vice versa
  _currentToProposedT = new double[_maxK];
  _proposedToCurrentT = new double[_maxK];
	

	
  _lambda = 0.5;
  _sigBeta = 0.05;
	
  _baseProbOfMove = 0.25;
  _baseProbOfDoNothing = 0.25;
	

}

Knots::~Knots()
{
  delete[] _currentK;

  //11/18/08 Chris Laws added following two lines to avoid memory leaks
  //shouldn't make much of a difference though because these are only called once jags is done
  for(unsigned int i = 0; i < _nchain; ++i)
    delete[] _currentT[i];
  delete[] _currentT;
  delete[] _proposedT;
	

  //11/18/08 Chris Laws added following two lines to avoid memory leaks
  //shouldn't make much of a difference though because these are only called once jags is done
  delete [] _currentToProposedT;
  delete [] _proposedToCurrentT;
	

}


void Knots::doNothing(unsigned int chain, RNG* const rng)
{
  
}

void Knots::birth(unsigned int chain, RNG* const rng)
{

  if(_currentK[chain] == _maxK)
    throw std::logic_error("cannot add knots if the number of knots is already at the max");
  
  //initialize the proposed values
  _proposedK = _currentK[chain];
  for(unsigned int i = 0; i < _proposedK; ++i)
    _proposedT[i] = _currentT[chain][i];

	
	
  //draw number of values to add
  _delta = rPoisT(_lambda, calcMaxDelta(_currentK[chain], Birth), rng);
	
	
  for(unsigned int i = 0; i < _delta; ++i)
    {
      //if there are currently no knots always propose from a uniform
      if(_currentK[chain] == 0)
        {
          _currentToProposedT[i] = rng->uniform();
        } 
      else
        {			
          //pick a position + 1 in the current knots to be the mean for the beta
          unsigned int T = static_cast<unsigned int>(rng->uniform() * (1 + _currentK[chain] - _minK));
			
          //if we picked position after all current knots then draw from a uniform
          if(T == _currentK[chain])
            _currentToProposedT[i] = rng->uniform();
          else
            _currentToProposedT[i] = rBeta(_currentT[chain][T], _sigBeta, rng);
        }

      addToProposedT(_currentToProposedT[i]);
    }

	
	
}


void Knots::death(unsigned int chain, RNG* const rng)
{
  if(_currentK[chain] == 0)
    throw std::logic_error("cannot remove knot if there are no knots");
  
  //initialize the proposed values
  _proposedK = _currentK[chain];
  for(unsigned int i = 0; i < _proposedK; ++i)
    _proposedT[i] = _currentT[chain][i];
	
  //draw number of values to remove
  _delta = rPoisT(_lambda, calcMaxDelta(_currentK[chain], Death), rng);
       
  for(unsigned int i = 0; i < _delta; ++i)
    {		
      //pick a position in the proposed knots to remove
      unsigned int T = static_cast<unsigned int>(rng->uniform() * (_proposedK - _minK));

      _proposedToCurrentT[i] = removeFromProposedT(T);
    }

}

void Knots::move(unsigned int chain, RNG* const rng)
{
  if(_currentK[chain] == 0)
    throw std::logic_error("cannot move knot if there are no knots");
  
  //initialize the proposed values
  _proposedK = _currentK[chain];
  for(unsigned int i = 0; i < _proposedK; ++i)
    _proposedT[i] = _currentT[chain][i];	
	
  //draw number of values to remove
  _delta = rPoisT(_lambda, calcMaxDelta(_currentK[chain], Death), rng);
	
  for(unsigned int i = 0; i < _delta; ++i)
    {
      //pick a position + 1 in the current knots to be the mean for the beta
      unsigned int T = static_cast<unsigned int>(rng->uniform() * (1 + _currentK[chain] - _minK));
			
      //if we picked position after all current knots then draw from a uniform
      if(T == _currentK[chain])
        _currentToProposedT[i] = rng->uniform();
      else
        _currentToProposedT[i] = rBeta(_currentT[chain][T], _sigBeta, rng);
	
    }
	
	
  for(unsigned int i = 0; i < _delta; ++i)
    {		
      //pick a position in the proposed knots to remove
      unsigned int T = static_cast<unsigned int>(rng->uniform() * (_proposedK - _minK));

      _proposedToCurrentT[i] = removeFromProposedT(T);
    }
	
  for(unsigned int i = 0; i < _delta; ++i)
    addToProposedT(_currentToProposedT[i]);


}

double Knots::acceptProbBalance(unsigned int const &chain, TypeOfUpdate type) const
{

  if(type == Nothing)
    throw std::logic_error("should never calculate \"acceptProbBalance\" for TypeOfUpdate equal to \"Nothing\"");

  double ans = 0;
	
  //K
  ans += std::log(1.0 / (_maxK - _minK)) - std::log(1.0 / (_maxK - _minK));
	
	
	
  //T
  for(unsigned int i = 0; i < _proposedK; ++i)
    ans += dbeta(_proposedT[i], _aPriorForT, _bPriorForT, 1);//std::log(1.0);
  for(unsigned int i = 0; i < _currentK[chain]; ++i)
    ans -= dbeta(_currentT[chain][i], _aPriorForT, _bPriorForT, 1);//std::log(1.0);
	
  //don't worry about delta because it cancels out
  if(type == Move)
    ans += std::log(calcProbOfMove(_proposedK)) - std::log(calcProbOfMove(_currentK[chain]));
	
  if(type == Birth)
    {
            
      //What is the prob of proposing backwards
      ans += std::log(calcProbOfDeath(_proposedK)) + dPoisT(_delta, _lambda, calcMaxDelta(_proposedK, Death));
            
            
      //What is the prob of proposing what we are proposing
      ans -= std::log(calcProbOfBirth(_currentK[chain])) + dPoisT(_delta, _lambda, calcMaxDelta(_currentK[chain], Birth));
            
    }
        
  if(type == Death)
    {
      //What is the prob of proposing backwards
      ans += std::log(calcProbOfBirth(_proposedK)) + dPoisT(_delta, _lambda, calcMaxDelta(_proposedK, Birth));
            
            
      //What is the prob of proposing what we are proposing
      ans -= std::log(calcProbOfDeath(_currentK[chain])) + dPoisT(_delta, _lambda, calcMaxDelta(_currentK[chain], Death));
            
    }
	
  if(type == Move)
    {
      //we don't have to worry about the ones which die off because that cancels out
      double tmp = 0.0;
            
      //calc prob of proposing current values from proposed values
      for(unsigned int i = 0; i < _delta; ++i)
        {
          tmp = 0.0;
          for(unsigned int j = 0; j < _proposedK; ++j)
            tmp += std::exp(dBeta(_proposedToCurrentT[i], _proposedT[j], _sigBeta)) / (_proposedK + 1);
          tmp += 1.0 / (_proposedK + 1);
		
          ans += std::log(tmp);
        }
            

      
      //calc prob of proposing porposed values
      for(unsigned int i = 0; i < _delta; ++i)
        {
          tmp = 0.0;
          for(unsigned int j = 0; j < _currentK[chain]; ++j)
            tmp += std::exp(dBeta(_currentToProposedT[i], _currentT[chain][j], _sigBeta)) / (_currentK[chain] + 1);
          tmp += 1.0 / (_currentK[chain] + 1);
		
          ans -= std::log(tmp);
        }
            
    }
	
  if(type == Birth)
    {
      double tmp = 0.0;
      
      //calc prob of proposing current values from proposed values
      
      for(unsigned int i = _proposedK; i >= _proposedK - _delta +1; --i)
        ans += std::log(1.0 / i);
            
            
      //calc prob of proposing porposed values
      for(unsigned int i = 0; i < _delta; ++i)
        {
          tmp = 0.0;
          for(unsigned int j = 0; j < _currentK[chain]; ++j)
            tmp += std::exp(dBeta(_currentToProposedT[i], _currentT[chain][j], _sigBeta)) / (_currentK[chain] + 1);
          tmp += 1.0 / (_currentK[chain] + 1);
		
          ans -= std::log(tmp);
        }
      
    }
        
  if(type == Death)
    {
            
      double tmp = 0.0;
            
      //calc prob of proposing current values from proposed values
      for(unsigned int i = 0; i < _delta; ++i)
        {
          tmp = 0.0;
          for(unsigned int j = 0; j < _proposedK; ++j)
            tmp += std::exp(dBeta(_proposedToCurrentT[i], _proposedT[j], _sigBeta)) / (_proposedK + 1);
          tmp += 1.0 / (_proposedK + 1);
			
          ans += std::log(tmp);
        }
		
		
      //calc prob of proposing porposed values
      for(unsigned int i = _currentK[chain]; i >= _currentK[chain] - _delta + 1; --i)
        ans -= std::log(1.0 / i); 
		
    }

  return ans;

}

Knots::TypeOfUpdate Knots::createProposal(RNG * const & r, unsigned int chain)
{
  
  unsigned int c = chain;
      
  double pDoNothing = calcProbOfDoNothing(_currentK[c]);
  double pMove = calcProbOfMove(_currentK[c]);
  
  //double pDeath = calcProbOfDeath(_currentK[c]);
  double pBirth = calcProbOfBirth(_currentK[c]);
  
  double rn = r->uniform();

  TypeOfUpdate type;
  if(rn < pDoNothing)
    {
      doNothing(c,r);
      type = Nothing;
      return type;
    } 
  else if (rn < pDoNothing + pMove)
    {
      move(c,r);
      type = Move;
      return type;
    }
  else if (rn < pDoNothing + pMove + pBirth)
    {
      birth(c,r);
      type = Birth;
      return type;
    }
  else 
    {
      death(c,r);
      type = Death;
      return type;
    }
  
}


unsigned int Knots::currentLength(unsigned int const &chain) const
{
  //return 1;
  return _currentK[chain];
}

double Knots::currentValue(unsigned int const &chain, unsigned int const & i) const
{
  //return 0.5 * (_maxT - _minT) + _minT;
  return _currentT[chain][i] * (_maxT - _minT) + _minT;
}

unsigned int Knots::proposedLength(unsigned int const &chain) const
{
  //return 1;
  return _proposedK;
}

double Knots::proposedValue(unsigned int const &chain, unsigned int const & i) const
{
  //return 0.5 * (_maxT - _minT) + _minT;
  return _proposedT[i] * (_maxT - _minT) + _minT;
}

double Knots::minT() const
{
  return _minT;
}

void Knots::acceptProposedValues(unsigned int const &chain)
{
  _currentK[chain] = _proposedK;

  for(unsigned int i = 0; i < _proposedK; ++i)
    _currentT[chain][i] = _proposedT[i];
}


double Knots::calcProbOfDoNothing(unsigned int K) const
{
  return _baseProbOfDoNothing;
}

double Knots::calcProbOfMove(unsigned int K) const
{
  if(K != 0)
    return _baseProbOfMove;
  else
    return 0;
}

double Knots::calcProbOfBirth(unsigned int K) const
{

  if(_minK == K)
    return 1.0 - calcProbOfDoNothing(K) - calcProbOfMove(K);

  if(K == _maxK)
    return 0.0;

  return 0.5 * (1.0 - calcProbOfDoNothing(K) - calcProbOfMove(K)) ;
}


double Knots::calcProbOfDeath(unsigned int K) const
{
  if(K == _minK)
    return 0.0;
	
  if(K == _maxK)
    return  1.0 - calcProbOfDoNothing(K) - calcProbOfMove(K);
	
  return 0.5 *(1.0 - calcProbOfDoNothing(K) - calcProbOfMove(K));
}
