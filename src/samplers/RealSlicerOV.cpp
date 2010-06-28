#include <graph/StochasticNode.h>
#include <graph/NodeError.h>
#include <distribution/Distribution.h>
#include <sampler/SampleMethod.h>
#include <RNG.h>
#include <util/nainf.h>

#include "RealSlicerOV.h"

#include <cmath>
#include <cfloat>

using std::vector;
using std::string;
using std::fabs;

#define MIN_ADAPT_OV 100
#define MIN_OV MIN_ADAPT_OV / 2



RealSlicerOV::RealSlicerOV(GraphView const * gv, unsigned int chain)
    : _gv(gv),
      _chain(chain),
      _probOfOverrelaxed(.9),
      _endpointAccuracy(5),
      _maxOV(100),
      _widthOV(1),
      _iterOV(0),
      _sumdiffOV(0),
      _adaptOV(true)
{

}

double RealSlicerOV::value() const
{
    return *(_gv->nodes()[0]->value(_chain));
}

void RealSlicerOV::setValue(double value) const
{
    _gv->setValue(&value, 1, _chain);
}

double RealSlicerOV::logFullConditional() const
{
    return _gv->logFullConditional(_chain);
}

bool 
RealSlicerOV::canSample(StochasticNode const *node)
{
    /*
      if (node->distribution()->isDiscreteValued() || node->length() != 1)
      return falses;
    
    
    if (node->distribution()->df(node) == 0)
	return false; 

    */

    if(node->distribution()->name() != "dchisqrOV" && 
       node->distribution()->name() != "dtOV" &&
       node->distribution()->name() != "dnormOV" &&
       node->distribution()->name() != "dunifOV")
	return false;
    
    return true;
}


void RealSlicerOV::update(RNG *rng)
{
    if(_iterOV < MIN_OV)
    {
	updateStep(rng);
    }
    else
    {
	if(rng->uniform() < _probOfOverrelaxed)
	    updateOverrelaxed(rng);
	else
	    updateStep(rng);
    }
}

void RealSlicerOV::getLimits(double *lower, double *upper) const
{
    _gv->nodes().front()->support(lower, upper, 1, _chain);
}


void RealSlicerOV::updateOverrelaxed(RNG *rng)
{
    // Test current value
    double g0 = logFullConditional();
    if (!jags_finite(g0)) {
	if (g0 > 0) {
	    return;
	}
	else {
	    throw NodeError(_gv->nodes()[0],
                            "Current value is inconsistent with data");
	}
    }

    // Generate auxiliary variable
    double z = g0 - rng->exponential();;

    // Generate random interval of width "_width" about current value
    const double xold = value();
    double L = xold - rng->uniform() * _widthOV; 
    double R = L + _widthOV;

    double lower = JAGS_NEGINF, upper = JAGS_POSINF;
    getLimits(&lower, &upper);

    // Stepping out 

    // Randomly set number of steps in left and right directions,
    // subject to the limit in the maximal size of the interval
    int j = static_cast<int>(rng->uniform() * _maxOV);
    int k = _maxOV - 1 - j;


    if (L < lower) {
	L = lower;
    }
    else {
	setValue(L);
	while (j-- > 0 &&logFullConditional() > z) {
	    L -= _widthOV;
	    if (L < lower) {
		L = lower;
		break;
	    }
	    setValue(L);
	}
    }

    if (R > upper) {
	R = upper;
    }
    else {
	setValue(R);
	while (k-- > 0 &&logFullConditional() > z) {
	    R += _widthOV;
	    if (R > upper) {
		R = upper;
		break;
	    }
	    setValue(R);
	}
    }

    /*Above portion is the same as updateStep*/
    /*Code below is the Overrelaxed portion*/

    double LBar = L;
    double RBar = R;
    double wBar = _widthOV;
    unsigned int aBar = _endpointAccuracy;

    if(R - L < 1.1 * _widthOV)
    {
	for(;;)
	{
	    double midPoint = (LBar + RBar) / 2;
	    setValue(midPoint);
	    if(aBar == 0 ||logFullConditional() > z)
		break;

	    if(xold > midPoint)
		LBar = midPoint;
	    else
		RBar = midPoint;

	    --aBar;
	    wBar /= 2.0;
	}
    }

    double LHat = LBar;
    double RHat = RBar;
    
    for(;aBar > 0; --aBar)
    {
	wBar /= 2.0;
	
	setValue(LHat + wBar);
	if(logFullConditional() <= z)
	    LHat += wBar;
	
	setValue(RHat - wBar);
	if(logFullConditional() <= z)
	    RHat -= wBar;	
    }

    double xnew = LHat + RHat - xold;

    setValue(xnew);
    if(xnew < LBar || xnew > RBar ||logFullConditional() <= z)
	setValue(xold);


    /*Code below portion is the same as updateStep*/
    if (_adaptOV && fabs(xnew - xold) > 0.0) {
	_sumdiffOV += _iterOV * fabs(xnew - xold);
	++_iterOV;
	if (_iterOV > MIN_ADAPT_OV) {
	    _widthOV = 2 * _sumdiffOV / _iterOV / (_iterOV - 1);  
	}
    }
}



void RealSlicerOV::updateStep(RNG *rng)
{
    // Test current value
    double g0 =logFullConditional();
    if (!jags_finite(g0)) {
	if (g0 > 0) {
	    return;
	}
	else {
	    throw NodeError(_gv->nodes()[0],
                            "Current value is inconsistent with data");
	}
    }

    // Generate auxiliary variable
    double z = g0 - rng->exponential();;

    // Generate random interval of width "_width" about current value
    double xold = value();
    double L = xold - rng->uniform() * _widthOV; 
    double R = L + _widthOV;

    double lower = JAGS_NEGINF, upper = JAGS_POSINF;
    getLimits(&lower, &upper);

    // Stepping out 

    // Randomly set number of steps in left and right directions,
    // subject to the limit in the maximal size of the interval
    int j = static_cast<int>(rng->uniform() * _maxOV);
    int k = _maxOV - 1 - j;


    if (L < lower) {
	L = lower;
    }
    else {
	setValue(L);
	while (j-- > 0 &&logFullConditional() > z) {
	    L -= _widthOV;
	    if (L < lower) {
		L = lower;
		break;
	    }
	    setValue(L);
	}
    }

    if (R > upper) {
	R = upper;
    }
    else {
	setValue(R);
	while (k-- > 0 &&logFullConditional() > z) {
	    R += _widthOV;
	    if (R > upper) {
		R = upper;
		break;
	    }
	    setValue(R);
	}
    }

    // Keep sampling from the interval until acceptance (the loop is
    // guaranteed to terminate).
    double xnew;
    for(;;) {
	xnew =  L + rng->uniform() * (R - L);
	setValue(xnew);
	double g =logFullConditional();
	if (g >= z - DBL_EPSILON) {
	    // Accept point
	    break;
	}
	else {
	    // shrink the interval
	    if (xnew < xold) {
		L = xnew;
	    }
	    else {
		R = xnew;
	    }
	}
    }

    if (_adaptOV && fabs(xnew - xold) > 0.0) {
	_sumdiffOV += _iterOV * fabs(xnew - xold);
	++_iterOV;
	if (_iterOV > MIN_ADAPT_OV) {
	    _widthOV = 2 * _sumdiffOV / _iterOV / (_iterOV - 1);  
	}
    }
}


string RealSlicerOV::name() const
{
    return "RealSlicerOV";
}


bool RealSlicerOV::adaptOff()
{
  _adaptOV = false;
  //FIXME We could try a bit harder than this
  return (_iterOV > MIN_ADAPT_OV );
}

bool RealSlicerOV::isAdaptive() const
{
    return true;
}




