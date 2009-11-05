#ifndef FUNCTION_H_
#define FUNCTION_H_

#include <string>
#include <vector>

class Range;

/**
 * Base class for functions.
 *
 * Function objects contain only constant data members and all member
 * functions are constant. Hence only one object needs to be
 * instantiated for each subclass.
 *
 * The FuncTab class provides a convenient way of storing Function
 * objects and referencing them by name.
 *
 * @see FuncTab
 * @short BUGS-language function or operator
 */
class Function
{
    const std::string _name;
    const unsigned int _npar;
public:
    /**
     * Constructor.
     * @param name name of the function as used in the BUGS language
     * @param npar number of parameters. If npar == 0 then the function
     * takes a variable number of parameters.
     */
    Function(std::string const &name, unsigned int npar);
    virtual ~Function ();
    /**
     * Returns the  BUGS language name of the function
     */
    std::string const &name () const;
    /**
     * Evaluates the function. 
     *
     * @param value array of doubles which contains the result of
     * the evaluation 
     * @param args Vector of arguments. 
     * @param lengths Respective lengths hof each element of args.
     * @param dims Respective dimensions of each element of args.
     */
    virtual void
	evaluate(double *value, 
		 std::vector<double const *> const &args,
		 std::vector<unsigned int> const &lengths,
		 std::vector<std::vector<unsigned int> > const &dims) const = 0;
    /**
     * Checks that a vector of parameters of length npar is consistent
     * with the function.
     */
    bool checkParameterLength(unsigned int npar) const;
    /**
     * Checks whether dimensions of the function parameters are correct.  
     *
     * @param dims Vector of length npar denoting the dimensions of
     * the parameters, with any redundant dimensions dropped.
     */
    virtual bool
	checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) 
	const = 0;
    /**
     * Checks whether the parameter values lie in the domain of the
     * function.  
     *
     * Since many functions do not need to check their parameter
     * values, a default method is provided which always returns
     * true.
     */
    virtual bool
	checkParameterValue(std::vector <double const *> const &args,
			    std::vector <unsigned int> const &lengths,
                            std::vector<std::vector<unsigned int> > const &dims) 
        const;
    /**
     * Calculates what the dimension of the return value should be,
     * based on the dimensions of the arguments. The dimension of the
     * return value cannot depend on the value of any of the arguments,
     * only their dimensions.
     *
     * This is a virtual function. The default implementation assumes
     * that the function is scalar-valued.
     *
     * @param dims Vector of Indices denoting the dimensions of the 
     * parameters. This vector must return true when passed to
     * checkParameterDim.
     */
    virtual std::vector<unsigned int> 
	dim(std::vector <std::vector<unsigned int> > const &dims) const;
    /**
     * Returns true if the function returns integer values. The default
     * implementation returns false, so any function that can return
     * integer values will need to overload this.
     *
     * @param mask Vector indicating whether parameters are discrete
     * or not. 
     *
     * @see SArray#isDiscreteValued
     */
    virtual bool isDiscreteValued(std::vector<bool> const &mask) const;
    /**
     * Tests whether the value of the function is linear, i.e. a
     * function of the form f(x) = A + B %*% x.  False negative
     * responses are permitted: i.e. the value false may be returned
     * when the value is, in fact, a linear function of the given
     * parameters, but false positives are not allowed. Since most
     * functions of interest are non-linear, the default method
     * returns false, so linear Functions must overload this.
     * 
     * @param mask boolean vector of length equal to the number of
     * parameters.  The mask indicates a subset of parameters (those with
     * value true) to be tested for linearity. The function may be 
     * non-linear for parameters with mask value false.
     *
     * @param isfixed boolean vector. This may be empty, in which case
     * it is ignored. A non-empty vector must have length equal to mask,
     * and denotes the parameters whose values are fixed. In this case
     * the test is more restrictive, and the function returns true only
     * if the coefficient B is fixed.
     *
     * @see DeterministicNode#isLinear
     */
    virtual bool isLinear(std::vector<bool> const &mask,
			  std::vector<bool> const &isfixed) const;
    /**
     * Tests whether the value of the function is a scale
     * transformation of the given scalar parameter, i.e. a trivial
     * linear function of the form f(x) = A + B %*% x, where A or B is
     * zero..  False negative responses are permitted.  Since most
     * functions of interest are not scale transformations, the
     * default method returns false.
     * 
     * @param index The index of the parameter to be tested. The function
     * may be a non-scale transformation of the other pararameters (e.g.
     * the function x/y is a scale transformation of x but not of y).
     *
     * @param isfixed boolean vector. This may be empty, in which case
     * it is ignored.  A non-empty vector must have length equal to
     * mask, and denotes the parameters whose values are fixed. In this
     * case the test is more restrictive, and the function returns true
     * only if the coeffiicent B is fixed.
     *
     * @see DeterministicNode#isScale
     */
    virtual bool isScale(unsigned int index,
			 std::vector<bool> const &isfixed) const;
    /**
     * Returns a BUGS-language expression representing the function call.
     * The default behaviour for a function named "foo" is to return
     * "foo(arg1,arg2)". Functions that are represented as prefix or infix
     * operators need to override this function.
     *
     * @param par Vector of parameter names for the function
     */
    virtual std::string deparse(std::vector<std::string> const &par) const;
};

#endif /* FUNCTION_H_ */
