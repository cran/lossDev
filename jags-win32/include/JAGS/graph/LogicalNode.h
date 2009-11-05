#ifndef LOGICAL_NODE_H_
#define LOGICAL_NODE_H_

#include <graph/DeterministicNode.h>
#include <vector>

class Function;

/**
 * @short Node defined by the BUGS-language operator <-
 *
 * A LogicalNode is a type of DeterministicNode that is defined by a
 * function, and a vector of parents as parameters. In the BUGS
 * language the function may be represented as an operator.
 */
class LogicalNode : public DeterministicNode {
    Function const * const _func;
    std::vector<std::vector<double const*> > _parameters;
    std::vector<unsigned int> _lengths;
    std::vector<std::vector<unsigned int> > _dims;
public:
    /**
     * A logical node is defined by a function (which may be an inline
     * operator in the BUGS language) and its parameters.
     */
    LogicalNode(Function const *func,
                std::vector<Node const*> const &parameters);
    ~LogicalNode();
    /** 
     * Calculates the value of the node based on the parameters. This
     * function calls Function#checkParameterValue and will throw an
     * exception if the result is false.
     * @exception NodeError
     */
    void deterministicSample(unsigned int chain);
    /**
     * @see Function#checkParameterValue.
     */
    bool checkParentValues(unsigned int chain) const;
    /**
     * A logical node preserves linearity if its function is linear.
     * @see Function#isLinear
     */
    bool isLinear(GraphMarks const &linear_marks, bool fixed) const;
    /**
     * A logical node preserves scale if its function is a scale function.
     * @see Function#isScale
     */
    bool isScale(GraphMarks const &scale_marks, bool fixed) const;
    std::string deparse(std::vector<std::string> const &) const;
};

#endif /* LOGICAL_NODE_H_ */
