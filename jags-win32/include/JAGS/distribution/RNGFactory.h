#ifndef RNG_FACTORY_H_
#define RNG_FACTORY_H_

#include <string>
#include <vector>

class RNG;

/**
 * @short Factory for RNG objects
 */
class RNGFactory
{
 public:
  /**
   * Destructor. An RNGFactory retains ownership of the RNG objects
   * it generates, and should delete them when the destructor is called.
   */
  virtual ~RNGFactory() {};
  /**
   * Returns a vector of newly allocated RNG objects.
   *
   * @param n Number of RNGs requested. If successful, the return value
   * will be a list of length n, and on exit the counter n will be
   * decremented to zero.  If a factory has the capacity to make only m < n
   * independent samplers, then the return value will be of length m
   * and on exit the counter will be decremented by m.
   */
  virtual std::vector<RNG *> makeRNGs(unsigned int &n) = 0;
  /**
   * Returns a newly allocated RNG object.  
   *
   * This function can be repeatedly called with the same name
   * argument. There is no guarantee that RNG objects created in this
   * way will generate independent streams.
   */
  virtual RNG * makeRNG(std::string const &name) = 0;
};

#endif /* RNG_FACTORY_H_ */
