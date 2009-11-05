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

#ifndef KNOTS_H_
#define KNOTS_H_
#include <JAGS/distribution/RNG.h>


/*
 * Values labeled as current are the current state.
 * We assume that at the begining of a call to update,
 * the current values contain the last accepted values.
*/
class Knots
{
 public:
	enum TypeOfUpdate {Nothing, Birth, Death, Move};

 protected:

	//min and max number of knots
	unsigned int _minK;
	unsigned int _maxK;
	
	//a and b for beta prior of Knots
	double _aPriorForT;
	double _bPriorForT;

	//min and max interval on which the Knots can fall
	//like in WinBUGS minT is also the position of Knot0 (which is not deemed to be part of _...T)
	double _minT;
	double _maxT;
	
	//number of knots
	unsigned int *_currentK;
	unsigned int _proposedK;
	
	
	//Position of Knots  these are stored on (0,1) and adjusted by _minT _maxT
	double **_currentT;
	double *_proposedT;
	
	//Values needed to go from _currentT to _proposedT and vice versa
	double *_currentToProposedT;
	double *_proposedToCurrentT;


        unsigned int _delta;
	double _lambda;
	double _sigBeta;
	
	double _baseProbOfMove;
	double _baseProbOfDoNothing;
        
	unsigned int _nchain;
        
	
private:
	
	double calcProbOfDoNothing(unsigned int K) const;
	double calcProbOfMove(unsigned int K) const;
	double calcProbOfBirth(unsigned int K) const;
	double calcProbOfDeath(unsigned int K) const;
	
	void doNothing(unsigned int chain,RNG *const rng);
	void birth(unsigned int chain, RNG *const rng);
	void death(unsigned int chain, RNG *const rng);
	void move(unsigned int chain, RNG *const rng);
	
	void addToProposedT(double const &T);
	double removeFromProposedT(unsigned int const &T);
	
	//calc max shift in K
	unsigned int calcMaxDelta(unsigned int const &K, TypeOfUpdate type) const;
	
	//draw from a Poisson Dist with mean lambda Truncated at T - 1 and return the answer + 1
	unsigned int rPoisT(double const &lambda, unsigned int const &T, RNG * const & rng) const;
	
	//return ll of getting x under rPoisT
	double dPoisT(unsigned int const &x, double const &lambda, unsigned int const &T) const;
	
	double rBeta(double const &mu, double const &sig, RNG * const & rng) const;
	double dBeta(double const &x, double const &mu, double const &sig) const;	
	
public:
	Knots(unsigned int nchain, double const * PriorForT, unsigned int const * KLimits, double const * TLimits);
	virtual ~Knots();
        double acceptProbBalance(unsigned int const &chain, TypeOfUpdate type) const;
        TypeOfUpdate createProposal(RNG * const & r, unsigned int chain);

        unsigned int currentLength(unsigned int const &chain) const;
        double currentValue(unsigned int const &chain, unsigned int const &i) const;

        unsigned int proposedLength(unsigned int const &chain) const;
        double proposedValue(unsigned int const &chain, unsigned int const &i) const;

        double minT() const;

        void acceptProposedValues(unsigned int const &chain);
        
};

#endif /*KNOTS_H_*/
