#ifndef DETERMINISTIC_NODE_H_
#define DETERMINISTIC_NODE_H_

#include <graph/Node.h>

/**
 * @short Base class for deterministic Node objects
 *
 * The value of a deterministic node is exactly by the values of its
 * parents.
 */
class DeterministicNode : public Node {
public:
    DeterministicNode(std::vector<unsigned int> const &dim,
		      std::vector<Node const *> const &parents);
    /**
     * Random samples from a Deterministic node are not random.
     * This function simply calculates the value of the node from its
     * parent values and leaves the RNG object untouched.
     */
    void randomSample(RNG*, unsigned int nchain);
    /**
     * Deterministic nodes are not random variables. This function
     * always returns false.
     */
    bool isRandomVariable() const;
};

#endif /* DETERMINISTIC_NODE_H_ */
