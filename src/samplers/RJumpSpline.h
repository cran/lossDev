/*
##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, 2010, 2011, 2012 National Council On Compensation Insurance Inc.,       ##
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

#ifndef RJUMPSPLINE_H_
#define RJUMPSPLINE_H_

#include <sampler/Sampler.h>
#include "Knots.h"
#include "ExpandableArray.h"
class DMNorm;

/*
 * Values labeled as current are the current state.
 * We assume that at the begining of a call to update,
 * the current values contain the last accepted values.
 */
class RJumpSpline : public Sampler
{

    GraphView * _mgv;
	
  //which node do does the sampler update
  StochasticNode const* _snode;

  Knots ** _knots;
  unsigned int _updatingKnotsI;
	
  //Coefficents
  //first value is an error term with zero mean and precision defined by the 8th parameter with must be NON-STOCHASTIC (_tauOfFirstIn3rdColumn)
  //next _TriDim-1 values are for the 3rd column of error terms after the first element
  //next _TriDim-1 values are for the 4th column of error terms (first value of this column is blank)
  //next value is intercept
  //next value is slope
  //next value is intercept
  //next value is slope
  // length of Coefficents is thus = 1 + 2 *(_TriDim-1) + 2 + 2 + number of knots + number of knots
  
  double **_currentBeta;
  double *_proposedBeta;
	
  //Prior for Betas and intercept
  double _tau;

  //Prior precision for 3rd and 4th column error terms;
  double _tauOf3rdColumn;
  double _tauOf4thColumn;

  //Prior precision for first value in the 3rd column NON-STOCHASTIC
  double _tauOfFirstIn3rdColumn;

  //autoregressive coefficient
  double _rho;
  double _etaRho;
  
  //values at which spline is to be evaluated
  //For now assume this is constant
  double **_x;
  //number of values at which spline is to be evaluated
  unsigned int *_n;

  unsigned int _TriDim;
	
	
  //posterior parameters for Beta under current and posterior cases
  double *_bPostCurrent;
  double *_APostCurrent;
	
  double *_bPostProposed;
  double *_APostProposed;
	
  //need this to calc ll of MNorm
  DMNorm* _dMNorm;
	
  unsigned int _nchain;
  unsigned int _numberOfSplines;

  ExpandableArray * _calPost_betas;
  ExpandableArray * _calPost_Acopy;
	
 private:
	
	
	
  //knotPositions are assumed to be in (0,1) and are scaled up to _minT _maxT
  void setSplineValue(bool const &current, unsigned int const &chain) const;
  
  unsigned int betaLength(bool const &current, unsigned int const& chain) const;

	
	
  //calc the ll of Z under current (or proposed) assume that double *_bPostCurrent, _APostCurrent, _bPostProposed, _APostProposed
  // have already be set and that betas are set to a likely value
  //this will change the value of the node
  double llZ(bool const &current, unsigned int chain) const;
	
  //these are taken out of JAGS/bug/MNormal
  //A call to either of these will change the value of the node and info stored in it will be lost
  void calBeta(double *betas, bool const &current, unsigned int const &chain);
  void calPost(bool const &current, unsigned int chain);
	
	
  void accept(unsigned int const &chain, unsigned int const &spline, Knots::TypeOfUpdate type, RNG * const & rng);
	

	
	
	
	
 public:
  RJumpSpline(GraphView *gv);
  virtual ~RJumpSpline();
  /**
   * Every sampler must update the vector of nodes and its immediate
   * deterministic descendants using the update function.
   *
   * @param rng vector of Pseudo-random number generator functions.
   */
  virtual void update(std::vector<RNG *> const &rng);
  /**
   * When a sampler is constructed, it may be in adaptive mode, which
   * allows it to adapt its behaviour for increased
   * efficiency. However, a sampler in adaptive mode may not converge
   * to the correct target distribution. This function turns off
   * adaptive mode, so that valid samples can be collected from the
   * sampler.
   */

  virtual void adaptOff();
  
  /** The adaptOff function may be called at any time. Premature ending
   * of adaptive mode may result in an extremely inefficient sampler.
   * The checkAdaptation function includes
   * an efficiency test to ensure that it has not been called
   * prematurely.  The return value is true if the efficiency test 
   * passes, and false otherwise.  Samplers that have no adaptive mode
   * should simply return true.
   */
  virtual bool checkAdaptation() const;
  /**
   * Indicates whether the sampler has an adaptive mode.
   */
  virtual bool isAdaptive() const;
  /**
   * Returns a name for the sampler which should describe the method
   * it uses to update the nodes.
   */
  virtual std::string name() const;
  /**
   * Static function that identifies the Marginal Stochastic Children
   * and the Immediate Deterministic Descendants of the given nodes
   * within the given graph.
   *
   * @param nodes Set of Nodes whose descendants are to be classified.
   *
   * @param graph Graph within which calculations are to take place.
   * Nodes outside of this graph will be ignored.
   *
   * @param stoch_nodes Vector which will contain the Marginal
   * Stochastic Children on exit.
   *
   * @param dtrm_nodes Vector which will contain the Immediate 
   * Deterministic Descendants, in forward sampling order, on exit.
   */
};

#endif /*RJUMPSPLINE_H_*/
