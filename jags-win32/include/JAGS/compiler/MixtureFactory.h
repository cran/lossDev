#ifndef MIXTURE_FACTORY_H_
#define MIXTURE_FACTORY_H_

#include <graph/Graph.h>
#include <graph/Node.h>

#include <vector>
#include <map>
#include <cfloat>

class NodeArray;
class MixtureNode;

/**
 * A MixMap is an STL map object.  The key represents an integer-valued
 * index, and the value is the corresponding Node from which a
 * MixtureNode will copy its values when its indices take that value.
 *
 * @see Mix
 */
typedef std::map<std::vector<int> , Node const *>  MixMap;

/**
 * A MixPair contains the information required to uniquely identify a
 * MixtureNode. The first element of the pair is a vector of Nodes
 * (which should be discrete-valued) representing the indices of the
 * MixtureNode. The second element is a MixMap which defines how each
 * possible value of the indices is mapped onto a parent Node.
 *
 * @see compMixPair ltmixpair
 */
typedef std::pair<std::vector <Node const*>, MixMap> MixPair;

/**
 * Comparison function for MixPair objects which is used to give them
 * a unique (but arbitrary) ordering. The compMixPair function returns
 * true if the first argument comes before the second in the ordering.
 *
 * @see ltmixpair
 */
bool compMixPair(MixPair const &, MixPair const &);

/**
 * @short STL function object for the map class using MixPair as a key
 *
 * @see MixtureFac
 */
struct ltmixpair
{
    bool operator()(MixPair const &arg1, MixPair const &arg2) const
    {
	return compMixPair(arg1, arg2);
    }
};

/**
 * @short Factory for MixtureNode objects
 * 
 * The purpose of a MixtureFactory is to avoid unnecessary duplication
 * of mixture nodes by having a container class and factory object
 * that will create and/or lookup mixture nodes.
 */
class MixtureFactory  
{ 
  std::map<MixPair, MixtureNode*, ltmixpair> _mixmap;
public:
  /**
   * Get a mixture node.  The results are cached, so if a request is
   * repeated, the same node will be returned.
   *
   * @param index Vector of discrete-valued Nodes that index the
   * possible parameters 
   *
   * @param parameters Vector of pairs. The first element of the pair
   * shows a possible value of the index vector and the second element
   * shows the corresponding parent from which the mixture node copies
   * its value.
   *
   * @param graph Graph to which newly allocated mixturenodes are added.
   *
   */
  MixtureNode *getMixtureNode(std::vector<Node const *> const &index,
			      MixMap const &parameters, Graph &graph);
};

#endif /* MIXTURE_FACTORY_H_ */
