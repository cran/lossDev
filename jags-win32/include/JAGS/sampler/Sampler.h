#ifndef SAMPLER_H_
#define SAMPLER_H_

#include <vector>
#include <set>
#include <string>

class StochasticNode;
class Node;
class Graph;
class RNG;

/**
 * @short Updates a set of stochastic nodes
 *
 * A sampler updates a set of stochastic nodes.  It is also
 * responsible for updating the immediate deterministic descendants of
 * those nodes (see below for a definition).  Sampling takes place in
 * the context of a Graph, which must contain the sampled nodes. Any
 * descendents of these nodes outside of the Graph are ignored when
 * updating.
 *
 * Some terminology:
 *
 * The "immediate deterministic descendants" of a set of stochastic
 * nodes S are the descendants of S in the graph where all stochastic
 * nodes except those in S have been removed.
 *
 * The "marginal stochastic children" of a set of stochastic nodes S
 * are the children of S in the graph where all deterministic nodes
 * have been marginalized out.
 *
 * A vector of nodes in an acyclic Graph is in "forward sampling
 * order" if node B always appears after node A when there is a path
 * from A to B. Note that forward sampling order is not uniquely
 * determined.
 */
class Sampler {
  unsigned int _length;
  std::vector<StochasticNode *> _nodes;
  std::vector<StochasticNode const *> _stoch_children;
  std::vector<Node*> _determ_children;
public:
  /**
   * Constructs a sampler for the given vector of nodes.  
   *
   * @param nodes Vector of Nodes to be sampled 
   *
   * @param graph  Graph within which sampling is to take place. It is
   * an error if this Graph does not contain all of the Nodes to be sampled.
   */
  Sampler(std::vector<StochasticNode *> const &nodes, Graph const &graph);
  virtual ~Sampler();
  /**
   * Returns the vector of sampled nodes
   */
  std::vector<StochasticNode *> const &nodes() const;
  /**
   * Sets the values of the sampled nodes.  The immediate
   * deterministic descendants are automatically updated.  This
   * function should be called by the update function.
   *
   * @param value Array of concatenated values to be applied to the sampled
   * nodes.
   *
   * @param length Length of the value array. This must be equal to the
   * sum of the  lengths of the sampled nodes.
   *
   * @param chain Number of the chain (starting from zero) to be modified.
   */
  void setValue(double const * value, unsigned int length, unsigned int chain);
  /**
   * Returns the total length of the sampled nodes
   */
  unsigned int length() const;
  /**
   * Returns the marginal stochastic children of the sampled nodes.
   */
  std::vector<StochasticNode const*> const &stochasticChildren() const;
  /**
   * Returns the immediate deterministic descendendants of the sampled
   * nodes, in forward sampling order
   */
  std::vector<Node*> const &deterministicChildren() const;
  /**
   * Calculates the log conditional density of the sampled nodes,
   * given all other nodes in the graph that was supplied to the
   * constructor, plus the parents of the nodes (which may be outside
   * the graph).  The log full conditional is calculated up to an
   * additive constant.
   *
   * @param chain Number of the chain (starting from zero) to query.
   */
  double logFullConditional(unsigned int chain) const;
  /**
   * Calculates the log prior density of the sampled nodes, i.e. the
   * density conditioned only on the parents.
   */
  double logPrior(unsigned int chain) const;
  /**
   * Calculates the log likelihood, which is added to the log prior
   * to give the log full conditional density
   */
  double logLikelihood(unsigned int chain) const;
  /**
   * Every sampler must update the vector of nodes and its immediate
   * deterministic descendants using the update function.
   *
   * @param rng vector of Pseudo-random number generator functions.
   */
  virtual void update(std::vector<RNG *> const &rng) = 0;
  /**
   * When a sampler is constructed, it may be in adaptive mode, which
   * allows it to adapt its behaviour for increased
   * efficiency. However, a sampler in adaptive mode may not converge
   * to the correct target distribution. This function turns off
   * adaptive mode, so that valid samples can be collected from the
   * sampler.
   *
   * The adaptOff function may be called at any time. Premature ending
   * of adaptive mode may result in an extremely inefficient sampler.
   * Therefore, any implementation of the adaptOff function must
   * include an efficiency test to ensure that it has not been called
   * prematurely.  The return value is true if the efficiency test 
   * passes, and false otherwise.  Samplers that have no adaptive mode
   * should simply return true.
   */
  virtual bool adaptOff() = 0;
  /**
   * Indicates whether the sampler has an adaptive mode.
   */
  virtual bool isAdaptive() const = 0;
  /**
   * Returns a name for the sampler which should describe the method
   * it uses to update the nodes.
   */
  virtual std::string name() const = 0;
  /**
   * Static function that identifies the Marginal Stochastic Children
   * and the Immediate Deterministic Descendants of the given nodes
   * within the given graph.
   *
   * @param nodes Set of Nodes whose descendants are to be classified.
   *
   * @param graph Graph within which calculations are to take place.
   * Nodes outside of this graph will be ignored.
   *
   * @param stoch_nodes Vector which will contain the Marginal
   * Stochastic Children on exit.
   *
   * @param dtrm_nodes Vector which will contain the Immediate 
   * Deterministic Descendants, in forward sampling order, on exit.
   */
  static void classifyChildren(std::vector<StochasticNode *> const &nodes,
			       Graph const &graph,
			       std::vector<StochasticNode const*> &stoch_nodes,
			       std::vector<Node*> &dtrm_nodes);
};

#endif /* SAMPLER_H_ */
