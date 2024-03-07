/**
 *    @file  basic_types_spec.hpp
 *   @brief  Specialisations for types defined in `basic_types.hpp`.
 *
 *  This header file defines macros for type specialisations particularly for
 *  generating ETIs.
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

#ifndef DIVERG_BASIC_TYPES_SPEC_HPP__
#define DIVERG_BASIC_TYPES_SPEC_HPP__

#define DIVERG_RANGEITERATOR_ETI_SPEC_AVAIL(ORDINAL)                       \
  extern template class RangeIterator< ORDINAL >;

#define DIVERG_RANGEITERATOR_ETI_SPEC_INST(ORDINAL)                        \
  template class RangeIterator< ORDINAL >;

#endif // DIVERG_BASIC_TYPES_SPEC_HPP__
