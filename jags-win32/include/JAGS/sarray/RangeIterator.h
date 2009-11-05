#ifndef RANGE_ITERATOR_H_
#define RANGE_ITERATOR_H_

#include <vector>

#include <sarray/Range.h>

/**
 * @short Mutable index that traverses a Range
 *
 * A RangeIterator is a numeric vector that is bound to be inside a given
 * Range.  It has operators to allow traversing the Range in row- or
 * column-major order.
 *
 * @see Range
 */
class RangeIterator : public std::vector<int> {
    const Range _range;
    unsigned int  _atend;
    //Forbid assignment
    RangeIterator &operator=(std::vector<int> const &);
public:
    /**
     * Constructor. The initial value of a RangeIterator is
     * the lower limit of the range argument.
     *
     * @param range. Range to traverse
     */
    RangeIterator(Range const &range);
    /**
     * Goes to the next index in column-major order, (i.e. moving the
     * left hand index fastest). If the RangeIterator at the upper
     * limit of the Range, then a call to nextLeft will move it to the
     * lower limit.
     * 
     * @return reference to self after incrementation
     * @see nextRight
     */
    RangeIterator &nextLeft();
    /**
     * Goes to the next index in row-major order (i.e. moving the right
     * hand index fastest) but otherwise behaves line nextLeft.
     *
     * @return reference to self after incrementation
     * @see nextLeft
     */
    RangeIterator &nextRight();
    /**
     * Returns a numeric counter of the number of times the
     * RangeIterator has gone past the upper bound of the Range (and
     * returned to the lower bound) in a call to nexLeft or nextRight.
     * The initial value is zero.
     */
    unsigned int atEnd() const;
};

#endif /* RANGE_ITERATOR_H_ */
