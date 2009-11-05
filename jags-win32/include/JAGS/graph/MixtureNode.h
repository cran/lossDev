#ifndef MIXTURE_NODE_H_
#define MIXTURE_NODE_H_

#include <graph/DeterministicNode.h>
#include <map>
#include <vector>

/**
 * @short Node for mixture models.
 *
 * A mixture node copies its value from one of several parents, based
 * on the current value of a vector of index nodes.
 *
 * In the BUGS language, mixture nodes are represented using nested
 * indexing. For example, in the deterministic relation
 *
 * <pre>y[i] <- x[ind[i]]</pre>
 *
 * y[i] is a mixture node if ind[i] is unobserved.  If the possible
 * values of ind[i] are 1...M, then the parents of y[i] are ind[i],
 * x[1], ... x[M].
 */
class MixtureNode : public DeterministicNode {
    std::map<std::vector<int>, Node const *> _map;
    unsigned int _Nindex;
public:
    /**
     * Constructs a MixtureNode. 
     *
     * @param index Vector of index nodes. These must be discrete-valued,
     *  scalar, and unobserved.
     *
     * @param mixmap STL map object which associates a possible value
     * of the index nodes with a single parent. Each possible index
     * value, denoted by a vector of integers, must have the correct
     * size (matching the size of the index parameter), and the
     * corresponding parents must all have the same dimension.
     */
    MixtureNode(std::vector<Node const *> const &index,
		std::map<std::vector<int>, Node const *> const &mixmap);
    ~MixtureNode();
    /**
     * Calculates the value of the mixture node by looking up the
     * parent node corresponding to the current value of the index
     * nodes (in the map supplied to the constructor). The mixture node
     * copies its value from that parent.
     */
    void deterministicSample(unsigned int chain);
    /**
     * Returns the number of index nodes.
     */
    unsigned int index_size() const;
    /**
     * A MixtureNode preserves linearity if none of the index nodes
     * are parameters. It is never a fixed linear function.
     */
    bool isLinear(GraphMarks const &linear_marks, bool fixed) const;
    /**
     * A MixtureNode is a scale transformation if none of the indices
     * are parameters. It is never a fixed scale transformation.
     */
    bool isScale(GraphMarks const &scale_marks, bool fixed) const;
    /**
     * This function always returns true. It ignores possible errors
     * in the index nodes.  This is because the deterministicSample
     * member function already checks that the index value is valid
     * and will throw an runtime error if it is not.  Repeating these
     * calculations in the checkParentValues function would therefore
     * be redundant.
     */
    bool checkParentValues(unsigned int chain) const;
    /**
     * Creates a name for the mixture node of the form.
     *
     * <pre>mixture(index=[p,q], parents=X[1,1,4]...X[2,3,4])</pre>
     *
     * Only the first and last parent nodes are listed to save space.
     *
     * Exceptionally, the name of a mixture node is not a reconstruction
     * of its BUGS-language definition.
     */
    std::string deparse(std::vector<std::string> const &parents) const;
};

bool isMixture(Node const *);
MixtureNode const * asMixture(Node const *);

#endif /* MIXTURE_NODE_H_ */
