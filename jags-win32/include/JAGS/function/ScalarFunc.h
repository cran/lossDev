#ifndef SCALAR_FUNC_H_
#define SCALAR_FUNC_H_

#include <function/Function.h>

/**
 * @short Scalar-valued Function with scalar arguments
 *
 * ScalarFunc is a convenience class for scalar-valued functions whose
 * parameters are also scalars. It provides simplified versions of some
 * of the virtual functions defined by the Function class. 
 *
 * Scalar functions are automatically vectorized so that, for example,
 * A + B are added element-wise if A and B are both conforming array
 * arguments. Scalar arguments are recycled so that, if c is a scalar,
 * A + c adds c to every element of A.
 */
class ScalarFunc : public Function
{
public:
    ScalarFunc(std::string const &name, unsigned int npar);
    /**
     * Implementation of the virtual function Function#evaluate, which
     * calls ScalarFunc#evaluateScalar.
     * 
     * This function automatically provides vectorization of the
     * scalar function defined by ScalarFunc#eval, applying it in turn
     * to elements of the argument vector "args". Scalar arguments are
     * recyled, as in the S language.
     * 
     * @param value Array of doubles which will contain the result on
     *              exit.
     * 
     * @param args Vector of arguments
     *
     * @param lengths Vector of argument lengths: the length of the array
     * of doubles pointed to by args[i] is lengths[i]. The lengths must
     * take one of two values: 1 or N where N is the size of the array
     * pointed to by the parameter "value".
     *
     * @dims Dimensions of arguments. These must be scalars or vectors
     * of length N, where N is the size of the 
     * 
     * not arrays.
     */
    void evaluate(double *value, 
		  std::vector <double const *> const &args,
		  std::vector <unsigned int> const &lengths,
                  std::vector<std::vector<unsigned int> > const &dims) const;
    /**
     * Scalar functions need to implement this simplified function,
     * instead of Function#evaluate.
     */
    virtual double evaluateScalar(std::vector <double const *> const &args) 
	const = 0;
    /**
     * Returns since the result of a scalar function should
     * always be scalar.
     */
    std::vector<unsigned int> 
	dim(std::vector <std::vector<unsigned int> > const &args) const;
    /**
     * Checks that all arguments are scalar, or that they are conforming
     * vector arguments.
     */
    bool checkParameterDim(std::vector<std::vector<unsigned int> > const &args)
	const;
    /**
     * Implementation of Function#checkParameterValue which calls the
     * simplified version ScalarFunc#checkScalarValue below. 
     */
    bool checkParameterValue(std::vector<double const *> const &args,
                             std::vector<unsigned int> const &lengths,
			     std::vector<std::vector<unsigned int> > const &dim)
	const;
    /**
     * Simplified version of checkParameterValue for functions that
     * take only scalar arguments.
     *
     * The default version always returns true.
     */
    virtual bool 
	checkScalarValue(std::vector<double const *> const &args) const;
};

#endif /* SCALAR_FUNC_H_ */
