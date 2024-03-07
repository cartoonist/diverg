/**
 *    @file  basic_types_impl.hpp
 *   @brief  Basic type definitions.
 *
 *  This header file defines some basic types.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Wed Mar 06, 2024  18:56
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2024, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_BASIC_TYPES_IMPL_HPP__
#define DIVERG_BASIC_TYPES_IMPL_HPP__

namespace diverg {
  /* https://stackoverflow.com/a/67642371/357257 */
  template< typename I >
  class RangeIterator
  {
    private:
      /* === DATA MEMBERS === */
      I i;
    public:
      typedef I difference_type;
      typedef I value_type;
      typedef I pointer;
      typedef I reference;
      typedef std::random_access_iterator_tag iterator_category;

      RangeIterator( I i ) : i( i ) { }

      bool operator==( RangeIterator<I>& other ) { return i == other.i; }
      I operator-( RangeIterator<I>& other ) { return i - other.i; }
      I operator++() { return i++; }
      I operator*() { return i; }
  };
}  /* --- end of namespace diverg --- */

#endif // DIVERG_BASIC_TYPES_IMPL_HPP__
