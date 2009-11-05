#ifndef MODEL_H_
#define MODEL_H_

#include <graph/Graph.h>

#include <vector>
#include <list>
#include <string>

class Sampler;
class Monitor;
class SamplerFactory;
class RNG;
class RNGFactory;
class MonitorFactory;

/**
 * @short Graphical model 
 *
 * The purpose of the model class is to collect together all the
 * elements necessary to run an MCMC sampler on a graphical model.
 */
class Model {
protected:
  std::vector<Sampler*> _samplers;
private:
  unsigned int _nchain;
  std::vector<RNG *> _rng;
  unsigned int _iteration;
  Graph _graph;
  std::set<Node*> _extra_nodes;
  std::vector<Node*> _sampled_extra;
  std::list<Monitor*> _monitors;
  std::list<Monitor*> _default_monitors;
  bool _is_initialized;
  bool _adapt;
  bool _data_gen;
  void initializeNodes(std::vector<Node*> const &sorted_nodes);
  void chooseRNGs();
  void chooseSamplers(std::vector<Node*> const &sorted_nodes);
  void setSampledExtra();
public:
  /**
   * @param nchain Number of parallel chains in the model.
   */
  Model(unsigned int nchain);
  virtual ~Model();
  /**
   * Returns the Graph associated with the model. This graph contains
   * all the nodes in the model
   */
  Graph &graph();
  /**
   * Initializes the model.  Initialization takes place in three steps.
   *
   * Firstly, random number generators are assigned to any chain that
   * doesn not already have an RNG.
   * 
   * Secondly, all nodes in the graph are initialized in forward
   * sampling order. 
   *
   * Finally, samplers are chosen for informative nodes in the graph.
   *
   * @param datagen Boolean flag indicating whether the model should
   * be considered a data generating model. If false, then
   * non-informative nodes will not be updated unless they are being
   * monitored.  This makes sampling more efficient by avoiding
   * redundant updates.  If true, then all nodes in the graph will be
   * updated in each iteration.
   *
   * @see Node#initialize, Model#rngFactories
   */
  void initialize(bool datagen);
  /** Returns true if the model has been initialized */
  bool isInitialized();
  /**
   * Updates the model by the given number of iterations. A
   * logic_error is thrown if the model is uninitialized.
   *
   * @param niter Number of iterations to run
   */
  void update(unsigned int niter);
  /**
   * Returns the current iteration number 
   */
  unsigned int iteration() const;
  /**
   * Adds a monitor to the model so that it will be updated at each
   * iteration.  This can only be done if Model#adaptOff has been
   * successfully called. Otherwise, a logic_error is thrown.
   */
  void addMonitor(Monitor *monitor);
  /**
   * Clears the monitor from the model, so that it will no longer
   * be updated. If the monitor has not previously been added to the
   * model, this function has no effect.
   */
  void removeMonitor(Monitor *monitor);
  /**
   * Returns the list of Monitors 
   */
  std::list<Monitor*> const &monitors() const;
  /**
   * Traverses the list of monitor factories requesting default
   * monitors of the given type. The function returns true after the
   * first monitor factory has added at least one node to the monitor
   * list. If none of the available monitor factories can create
   * default monitors of the given type, the return value is false.
   *
   * Monitors created by a call to setDefaultMonitors are owned by
   * the Model.
   *
   * @see MonitorFactory#addDefaultMonitors
   */
  bool setDefaultMonitors(std::string const &type, unsigned int thin);
  /**
   * Removes all Monitors of the given type created by a previous call
   * to setDefaultMonitors, and deletes them.
   */
  void clearDefaultMonitors(std::string const &type);
  /**
   * After the model is initialized, extra uninformative nodes may be
   * added to the graph (For example, a node may be added that
   * calculates the deviance of the model).  The model takes
   * responsibility for updating the extra node.
   *
   * The extra node cannot be observed, it must not already be in the
   * model graph, it may not have any children, and all of its parents
   * must be in the graph.
   */
  void addExtraNode(Node *node);
  /**
   * Access the list of sampler factories, which is common to all
   * models. This is used during initialization to choose samplers.
   *
   * @seealso Model#chooseSamplers
   */
  static std::list<SamplerFactory const *> &samplerFactories();
  /**
   * Access the list of RNG factories, which is common to all models
   */
  static std::list<RNGFactory *> &rngFactories();
  /**
   * Access the list of monitor factories, which is commmon to all models
   */
  static std::list<MonitorFactory *> &monitorFactories();
  /**
   * Returns the number of chains in the model
   */
  unsigned int nchain() const;
  /**
   * Returns the RNG object associated with the given chain. If no RNG
   * has been assigned, then a NULL pointer is returned.
   */
  RNG *rng(unsigned int nchain) const;
  /**
   * Assigns a new RNG object to the given chain. The list of
   * RNGFactory objects is traversed and each factory is requested to
   * generate a new RNG object of the given name.
   *
   * @return success indicator.
   */
  bool setRNG(std::string const &name, unsigned int chain);
  /**
   * Assigns an existing RNG object to the given chain
   *
   * @return success indicator
   */
  bool setRNG(RNG *rng, unsigned int chain);
  /**
   * Turns off adaptive phase of all samplers.
   *
   * @return true if all samplers passed the efficiency
   * test. Otherwise false.
   *
   * @see Sampler#adaptOff
   */
  bool adaptOff();
  /**
   * Indicates whether the model is in adaptive mode (before the
   * adaptOff function has been called).
   */
  bool isAdapting() const;
};

#endif /* MODEL_H_ */
