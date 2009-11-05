#ifndef CONSTANT_FACTORY_H_
#define CONSTANT_FACTORY_H_

#include <map>
#include <compiler/NodeFactory.h>

class ConstantNode;
class Graph;

/**
 * @short STL function object for the map class using double as a key
 */
struct ltdouble
{
  bool operator()(double arg1, double arg2) const
  {
    return lt(arg1, arg2);
  }
};

typedef std::pair<std::vector<unsigned int>, std::vector<double> > constpair;

/**
 * @short Factory for ConstantNode objects
 *
 * The purpose of a ConstantFactory is to avoid unnecessary
 * duplication of constant nodes by having a container class and
 * factory for them that will create and/or lookup constant
 * nodes based on their value.
 */
class ConstantFactory 
{ 
    unsigned int _nchain;
    std::map<double, ConstantNode*, ltdouble> _constmap;
    std::map<constpair, ConstantNode*> _mv_constmap;
public:
    ConstantFactory(unsigned int nchain);
    /**
     * Get a constant node with a given value.  The results are cached,
     * so if a request is repeated, the same node will be returned.
     * If a node is newly allocated, it is inserted into the given graph.
     */
    ConstantNode *getConstantNode(double value, Graph &graph);
    /**
     * Get a multivariate constant node. The results are cached
     * so that if a request is repeated, the same node will be returned.
     * If a node is newly allocated, it is inserted into the given graph.
     */
    ConstantNode *getConstantNode(std::vector<unsigned int> const &dim,
				  std::vector<double> const &value,
				  Graph &graph);
};

#endif /* CONSTANT_FACTORY_H_ */
