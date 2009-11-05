#ifndef MONITOR_FACTORY_H_
#define MONITOR_FACTORY_H_

#include <string>
#include <vector>

class Monitor;
class Node;
class Model;

/**
 * @short Factory for Monitor objects
 */
class MonitorFactory {
public:
    virtual ~MonitorFactory();
    /**
     * Creates a monitor of the given type. If the monitor cannot
     * be constructed, a null pointer is returned
     */
    virtual Monitor *
	getMonitor(Node const *node, Model *model, unsigned int start,
		   unsigned int thin, std::string const &type) = 0;
    /**
     * Returns a vector of default nodes. These are nodes for which
     * monitors are typically created manually by the user.
     */
    virtual std::vector<Node const*> defaultNodes(Model *model,
						  std::string const &type) 
	const = 0;
};

#endif /* MONITOR_FACTORY_H_ */
