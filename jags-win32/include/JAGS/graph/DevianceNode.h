#ifndef DEVIANCE_NODE_H_
#define DEVIANCE_NODE_H_

#include <graph/DeterministicNode.h>
#include <graph/StochasticNode.h>
#include <vector>
#include <set>
#include <stdexcept>

class Model;
class StochasticNode;

/**
 * @short Special node to calculate deviance statistic
 */
class DevianceNode : public DeterministicNode
{
    std::vector<StochasticNode const*> _parameters;
public:
    /**
     * A DevianceNode is a scalar node whose parents are all
     * stochastic nodes. The value of a Deviance node is minus twice
     * the log density of the parents.
     */
    DevianceNode(std::set<StochasticNode const*> const &parameters);
    /**
     * Calculate the deviance.
     */
    void deterministicSample(unsigned int chain);
    /**
     * A deviance node never preserves linearity.
     */
    bool isLinear(GraphMarks const &linear_marks, bool fixed) const;
    /**
     * A deviance node is never a scale transformation of its parameters.
     */
    bool isScale(GraphMarks const &scale_marks, bool fixed) const;
    /**
     * A deviance node always has valid parent values
     */
    bool checkParentValues(unsigned int chain) const;
    std::string deparse(std::vector<std::string> const &parents) const;
};

#endif
