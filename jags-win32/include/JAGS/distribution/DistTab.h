#ifndef DIST_TAB_H_
#define DIST_TAB_H_

#include <list>
#include <string>

class Distribution;

/** 
 * @short Look-up table for Distribution objects
 *
 * Since all member functions of the Distribution class are constant,
 * only one instance of a Distribution object is required. The DistTab
 * class provides a means of storing Distributions and looking them up
 * by name.
 *
 * @see FuncTab
 */
class DistTab {
    std::list<Distribution const *> _dist_list;
    std::list<Distribution const *> _masked_list;
public:
  /**
   * Inserts a distribution into the table.  Adding a new distribution
   * with the same name as an existing distribution in the table has
   * the effect of masking (but not erasing) the existing entry.
   *
   * @param dist Pointer to the distribution to insert. The pointer
   * must be valid for the lifetime of the DistTab, unless the
   * distribution is removed with a call to erase.
   */
  void insert(Distribution const *dist);
  /**
   * Finds a distribution by name. If two distributions with the
   * same name have been inserted into the table, the most recently
   * inserted one is returned.
   *
   * @return a pointer to the distribution or 0 if it was not found.
   */
  Distribution const *find(std::string const &name) const;
  /**
   * Removes a distribution from the table. A DistTab does no memory
   * management and will not delete the erased distribution.  This
   * must be done by the owner of the distribution. Conversely only
   * the owner can call the erase function since it requires a non-constant
   * pointer.
   */
  void erase(Distribution *dist);
  /**
   * Returns the list of Distributions accessible by the find function.
   */
  std::list<Distribution const *> const &distributions() const;
  /**
   * Returns the list of distributions that are currently masked by a
   * different distribution with the same name.
   */
  std::list<Distribution const *> const &masked() const;
};

#endif /* DIST_TAB_H_ */
