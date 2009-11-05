#ifndef NODE_H_
#define NODE_H_

#include <set>
#include <vector>
#include <string>

class RNG;
class Graph;
class GraphMarks;

/**
 * @short Node in a directed acyclic graph representingn a Bayesian model 
 *
 * Nodes are reference managed and will delete themselves when the
 * reference count reaches zero.  Referencing and dereferencing
 * takes place automatically when a node is inserted into/removed
 * from a graph, but can also be done by calling the ref and
 * unref member functions.
 *
 * Nodes should be dynamically allocated, and should be inserted
 * into a Graph as soon as they are created.  
 *
 * @see Graph
 */
class Node {

  double **_dataByChain;


    std::vector<Node const *> _parents;
    std::set<Node*> *_children;
    unsigned int _ref;
    bool _isobserved;
    bool _isdiscrete;

    /* Forbid copying of Node objects */
    Node(Node const &orig);
    Node &operator=(Node const &rhs);
protected:
    const std::vector<unsigned int> _dim;
    const unsigned int _length;
    const unsigned int _nchain;
    double *_data;
  
public:
    /**
     * Constucts a Node with no parents.
     * @param dim Dimension of new Node.
     * @param nchain Number of chains that the Node will contain data on.
     */
    Node(std::vector<unsigned int> const &dim, unsigned int nchain);
    /**
     * Constructs a node with parents.  Each parent must contain data on
     * the same number of chains. Subclasses of Node may give specific
     * meaning to the ordering of the parents.
     *
     * @param dim Dimension of new Node.
     *
     * @param parents vector of parent nodes. A node may not be its own
     * parent.
     */
    Node(std::vector<unsigned int> const &dim, 
	 std::vector<Node const *> const &parents);
    /**
     * Destructor. 
     */
    virtual ~Node();
    /**
     * Increments the reference count.
     */
    void ref();
    /**
     * Decrements reference count. The node deletes itself when the
     * reference count is zero.
     */
    void unref();
    /**
     * Shows the current reference count.
     */
    unsigned int refCount() const;
    /**
     * Number of chains.
     */ 
    unsigned int nchain() const;
    /**
     * Vector of parents.
     */
    std::vector<Node const *> const &parents() const;
    /**
     * Set of children.
     *
     * Note that if we have write access to a Node, then this function
     * gives write access to it's children.  This is a necessity:
     * if we modify the value of the current Node, then we may need to
     * update it's children to keep consistency of the model. Conversely,
     * if we do not have write access to the Node (e.g. we have a constant
     * pointer or reference) then this function cannot be used to obtain
     * access to it's children.
     */
    std::set<Node*> const *children();
    /**
     * Draws a random sample from the node's prior distribution.
     * @param rng Pointer to random number generator
     * @param chain Number of chain from which to draw sample
     */
    virtual void randomSample(RNG *rng, unsigned int chain) = 0;
    /**
     * Calculates a value for the node based on its parents' values.
     * @param chain Number of chain from which to draw sample
     */
    virtual void deterministicSample(unsigned int chain) = 0;
    /**
     * Checks whether the parents of the Node have valid values.
     */
    virtual bool checkParentValues(unsigned int chain) const = 0;
    /**
     * Initializes the node for the given chain. The value array of a
     * newly constructed Node consists of missing values (denoted by
     * the special value JAGS_NA).  This function sets the value of
     * the node by forward sampling from its parents.  If the Node has
     * previously had its value set, the function will do nothing and
     * return the value true.  Initialization will fail if any of the
     * parent nodes is uninitialized, and in this case the return
     * value is false.
     *
     * @param rng Random number generator 
     * 
     * @param chain Index number of chain to initialize.
     *
     * @returns a logical value indicating success
     */
    bool initialize(RNG *rng, unsigned int chain);
    /**
     * Initializes a node, in all chains, if it is not a random
     * variable and if all of its parents are observed. In this case,
     * the value of the node is also fixed. Otherwise the function
     * has no effect.
     *
     * @see initialize
     */
    void initializeData();
    /**
     * Returns the BUGS-language representation of the node, based on the 
     * names of its parents
     *
     * @param parents Vector of names of parent nodes
     */
    virtual std::string deparse(std::vector<std::string> const &parents) const = 0;
    /**
     * Returns true if the node represents a random variable.
     */
    virtual bool isRandomVariable() const = 0;
    /**
     * Sets the value of the node, in all chains, and marks the node
     * as observed. This function can only be called once for a given
     * Node.
     * 
     * @param value vector of values
     */
    void setObserved(std::vector<double> const &value);
    /**
     * Indicates whether the node is observed. 
     */
    bool isObserved() const;
    /**
     * Sets the value of the node for a given chain
     * @param value Array of values to be assigned
     * @param length Length of the value array
     * @param chain number of chain (starting from zero) to modify
     *
     * @see SArray#setValue
     */
    void setValue(double const *value, unsigned int length, unsigned int chain);
    /**
     * Indicates whether a node is discrete-valued or not.
     * @see SArray#isDiscreteValued
     */
    bool isDiscreteValued() const;
    /**
     * Permanently sets the node to be discrete-valued for all chains.
     * @see SArray#isDiscreteValued
     */
    void setDiscreteValued();
    /**
     * Returns a pointer to the start of the array of values for 
     * the given chain.
     */
    double const *value(unsigned int chain) const;
    /**
     * Returns the length of the value array
     */
    unsigned int length() const;
    /**
     * Returns the dimensions of the Node
     */
    std::vector<unsigned int> const &dim() const;
    /**
     * Tests whether the value of the node is a linear function of a
     * set of ancestor nodes X = (X1, ... Nn), i.e. whether the value
     * of the node can be expressed as A + B %*% X1 + B2 %*% X2 + ...
     * + Bn %*% Xn. Preservation of linearity is a criterion used by
     * some Samplers to determine if they can act on a set of
     * stochastic nodes.
     *
     * False negative responses are permitted: i.e. the value false may
     * be returned when the node is, in fact, a linear function, but
     * false positives are not allowed.
     * 
     * The test for linearity takes place in a graph, which is implicitly
     * defined by the linear_marks parameter (and may be accessed directly
     * using the GraphMarks#graph member function on linear_marks). Only
     * paths that are contained entirely inside the graph are considered.
     * The graph is assumed to be acyclic. Since the Node#isLinear function
     * is designed to be called iteratively on a sequence of nodes, it
     * relies on a GraphMarks object for book-keeping.
     *
     * @param linear_marks A GraphMarks object in which all ancestors
     * of the current node that are also descendants of X have been
     * marked.  The mark values are MARK_TRUE if the node is a
     * (possibly fixed) linear function of X, and MARK_FALSE
     * otherwise.  Ancestors of the current node that are not
     * descendants of X should remain unmarked.
     *
     * @param fixed Logical flag. When true, the test is more
     * stringent and the isLinear function returns the value true only
     * if the current node is a fixed linear function of X, i.e. the
     * coefficients B1, ... Bn are constant (but not necessarily
     * A). In this case, the linear_marks parameter must also conform
     * to the more stringent conditions: ancestors of the current node
     * should be marked with MARK_TRUE only if they are fixed linear
     * functions of X. They should be marked with MARK_FALSE if they
     * are non-linear or non-fixed linear functions of X.
     */
    virtual bool 
	isLinear(GraphMarks const &linear_marks, bool fixed) const = 0;
    /**
     * Tests whether the value of the node is a scale function of the
     * ancestor node X. A scale function is a trivial linear function
     * of the form A + B %*% X where either A or B is zero (or both).
     * Preservation of scale is used by some Samplers to determine
     * if they can act on a set of stochastic nodes.
     * 
     * The isScale function works the same way as the isLinear function.
     *
     * @param scale_marks. GraphMarks object in which all ancestors
     * of the current node that are also descendants of X are marked
     * with MARK_TRUE or MARK_FALSE, depending on whether they are
     * (possibly fixed) scale functions of X or not.
     *
     * @param fixed When true, the test is more stringent and returns
     * the value true only if the function is a scale transformation
     * with fixed coefficient B. In this case, the scale_marks parameter
     * must conform to the more stringent conditions: only fixed scale
     * functions of X may be marked with MARK_TRUE
     */
    virtual bool 
	isScale(GraphMarks const &scale_marks, bool fixed) const = 0;
};

/**
 * Calculate the number of chains of parameter vector. Returns 0
 * if the parameters are inconsistent
 */
unsigned int countChains(std::vector<Node const *> const &parameters);

#endif /* NODE_H_ */
