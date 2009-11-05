#ifndef NODE_FACTORY_H_
#define NODE_FACTORY_H_

#include <graph/Graph.h>

#include <vector>
#include <cfloat>

inline bool lt(double value1, double value2)
{
    return value1 < value2 - 16 * DBL_EPSILON;
}

bool lt(std::vector<double> const &value1, std::vector<double> const &value2);
bool lt(double const *value1, double const *value2, unsigned int length);
bool lt(Node const *node1, Node const *node2);
bool lt(std::vector<Node const *> const &par1, 
	std::vector<Node const *> const &par2);

#endif /* NODE_FACTORY_H_ */
