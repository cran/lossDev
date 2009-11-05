#ifndef BUGS_MODEL_H_
#define BUGS_MODEL_H_

#include <vector>
#include <map>
#include <utility>
#include <string>

#include <model/Model.h>
#include <model/SymTab.h>

/**
 * @short Model with symbol table 
 *
 * A BUGS model is a subclass of Model that contains a symbol table
 * which is used to store certain nodes in arrays.  The array format
 * gives a convenient way of looking up nodes by name.
 */
class BUGSModel : public Model
{
    SymTab _symtab;
    std::map<Node const*, std::pair<std::string, Range> > _node_map;
    std::list<Monitor*> _bugs_monitors;
public:
    BUGSModel(unsigned int nchain);
    ~BUGSModel();
    /**
     * Returns the symbol table of the BUGSModel.
     */
    SymTab &symtab();
    /**
     * Adds a DevianceNode to the model graph.  The parents of this node
     * are all observed stochastic nodes in the graph.  A new NodeArray
     * is added to the symbol table with name "deviance" and the new
     * DevianceNode is inserted into it. However, if a NodeArray with
     * name "deviance" already exists in the symbol table then this
     * function does nothing.
     */
    void addDevianceNode();
    /**
     * Returns a node given a variable name and a range of indices. The
     * name is used to look up a NodeArray in the symbol table, and 
     * a subset of the Nodearray corresponding to the given index range
     * is taken. 
     *
     * If the requested Node cannot be returned, a NULL pointer is 
     * returned and an informative error message is given.
     *
     * @param name String giving the BUGS-language name of a variable
     *
     * @param range Range representing the lower and upper limits of the
     * subset. A NULL range corresponds to the whole range of the
     * NodeArray.
     *
     * @param message If the requested Node does not exist, an
     * informative error message is copied to this string.
     */
    Node *getNode(std::string const &name, Range const &range,
		  std::string &message);
    /**
     * Writes out stored samples in CODA format.
     *
     * @param nodes Vector of nodes to write out. Each node is described
     * by a pair consisting of a variable name and a range of indices.
     * If a node is not being monitored, then it is ignored.
     *
     * @param index ofstream for the index file. The index file consists
     * of three entries per line: the node name, the start line number,
     * and the end line number. Start and end line numbers refer to the
     * lines of the output file(s). A multivariate node has separate
     * entries for each element, e.g. A[1,1], A[2,1], etc.
     * 
     * @param output vector of ofstreams for the output files. There
     * must be one output file for each chain in the model or a
     * logic_error is thrown. An output file consists of two entries per
     * line: the iteration number and the node value.
     *
     * @exception logic_error
     */
    void coda(std::vector<std::pair<std::string,Range> > const &nodes, 
	      std::ofstream &index, std::vector<std::ofstream*> const &output,
              std::string &warn);
    void coda(std::ofstream &index, std::vector<std::ofstream*> const &output,
              std::string &warn);
    /**
     * Sets the state of the RNG, and the values of the unobserved
     * stochastic nodes in the model, for a given chain.
     *
     * @param param_table STL map, in which each entry relates a
     * variable name to an SArray.  If the name is ".RNG.state" or
     * ".RNG.seed", then the SArray is used to set the state of the RNG.
     * Otherwise the SArray value is used to set the unobserved
     * stochastic nodes in the model. Each SArray must have the same
     * dimensions as the NodeArray in the symbol table with the
     * corresponding name. Elements of the SArray must be set to the
     * missing value, JAGS_NA, unless they correspond to an element of
     * an unobserved StochasticNode.
     *
     * @param chain Number of chain (starting from zero) for which
     * parameter values should be set.
     *
     * @see RNG#init RNG#setState
     * @exception runtime_error
     */
    void setParameters(std::map<std::string, SArray> const &param_table,
		       unsigned int chain);
    /**
     * Creates a new Monitor. The BUGSModel is responsible for the
     * memory management of any monitor created this way. It is not
     * possible to create two monitors with the same name, range and
     * type.
     *
     * @param name Name of the node array
     *
     * @param range Subset of indices of the node array defining hte
     * node to be monitored.
     * 
     * @param thin Thinning interval for monitor
     *
     * @param type Type of monitor to create
     *
     * @return True if the monitor was created.  
     */
    bool setMonitor(std::string const &name, Range const &range,
		    unsigned int thin, std::string const &type);
    /**
     * Deletes a Monitor that has been previously created with a call
     * to setMonitor.
     *
     * @return True if the monitor was deleted.
     */
    bool deleteMonitor(std::string const &name, Range const &range,
		       std::string const &type);
    /**
     * Writes the names of the samplers, and the corresponding 
     * sampled nodes vectors to the given vector.
     *
     * @param sampler_names vector that is modified during the call On
     * exit it will contain an element for each Sampler in the model.
     * Each element is a vector of strings: the first string is the
     * name of the sampler, and the remaining strings are the names of
     * the nodes sampled by that Sampler.
     */
    void samplerNames(std::vector<std::vector<std::string> > &sampler_names) 
	const;

};

#endif /* BUGS_MODEL_H_ */
