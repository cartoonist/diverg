/**
 *    @file  bits.hpp
 *   @brief  Interface functions for bitwise operations.
 *
 *  This header file contains general utility and helper functions for bitwise
 *  operations.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Tue Mar 12, 2024  16:01
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2024, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_BITS_HPP__
#define DIVERG_BITS_HPP__

#include <stddef.h>
#include <stdint.h> // for uint64_t uint32_t declaration


namespace diverg {
  /**
   *   @brief  A helper class for bitwise tricks on 64 bit words.
   *
   *   NOTE: This class is adopted from `sdsl::bits.hpp` header file.
   *
   *   ---
   *   Notes from original implementation:
   *   bits is a helper class for bitwise tricks and techniques.
   *   For the basic tricks and techiques we refer to Donald E. Knuth's
   *   "The Art of Computer Programming", Volume 4A, Chapter 7.1.3 and
   *   the informative website of Sean E. Anderson about the topic:
   *   http://www-graphics.stanford.edu/~seander/bithacks.html .
   *
   *   We have added new functions like: cnt11 and sel11.
   *
   *   All members of this class are static variables or methods.
   *   This class cannot be instantiated.
   *
   *   ---
   *
   *   @author Simon Gog
   */
  template< typename T = void >
  struct bits_impl
  {
    bits_impl() = delete;

    /**
     *   @brief This constant represents a de Bruijn sequence B(k,n) for k=2 and n=6.
     *
     *   Details for de Bruijn sequences see:
     *   http://en.wikipedia.org/wiki/De_bruijn_sequence
     *   deBruijn64 is used in combination with the array lt_deBruijn_to_idx.
     */
    static constexpr uint64_t deBruijn64{0x0218A392CD3D5DBFULL};

    /**
     *   @brief This table maps a 6-bit subsequence S[idx...idx+5] of constant deBruijn64 to idx.
     *   @sa deBruijn64
     */
    static constexpr uint32_t lt_deBruijn_to_idx[64] = {0,  1,  2,  7,  3,  13, 8,  19, 4,  25, 14, 28, 9,  34, 20, 40,
                                                        5,  17, 26, 38, 15, 46, 29, 48, 10, 31, 35, 54, 21, 50, 41, 57,
                                                        63, 6,  12, 18, 24, 27, 33, 39, 16, 37, 45, 47, 30, 53, 49, 56,
                                                        62, 11, 23, 32, 36, 44, 52, 55, 61, 22, 43, 51, 60, 42, 59, 58};

    /** Lookup table for byte popcounts. */
    static constexpr uint8_t lt_cnt[256] = {
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2,
        3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3,
        3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
        6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4,
        3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4,
        5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6,
        6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

    /** Lookup table for most significant set bit in a byte. */
    static constexpr uint32_t lt_hi[256] = {
        0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

    /** lo_set[i] is a 64-bit word with the i least significant bits set and the high bits not set. */
    /*! lo_set[0] = 0ULL, lo_set[1]=1ULL, lo_set[2]=3ULL...
     */
    static constexpr uint64_t lo_set[65] = {
        0x0000000000000000ULL, 0x0000000000000001ULL, 0x0000000000000003ULL, 0x0000000000000007ULL,
        0x000000000000000FULL, 0x000000000000001FULL, 0x000000000000003FULL, 0x000000000000007FULL,
        0x00000000000000FFULL, 0x00000000000001FFULL, 0x00000000000003FFULL, 0x00000000000007FFULL,
        0x0000000000000FFFULL, 0x0000000000001FFFULL, 0x0000000000003FFFULL, 0x0000000000007FFFULL,
        0x000000000000FFFFULL, 0x000000000001FFFFULL, 0x000000000003FFFFULL, 0x000000000007FFFFULL,
        0x00000000000FFFFFULL, 0x00000000001FFFFFULL, 0x00000000003FFFFFULL, 0x00000000007FFFFFULL,
        0x0000000000FFFFFFULL, 0x0000000001FFFFFFULL, 0x0000000003FFFFFFULL, 0x0000000007FFFFFFULL,
        0x000000000FFFFFFFULL, 0x000000001FFFFFFFULL, 0x000000003FFFFFFFULL, 0x000000007FFFFFFFULL,
        0x00000000FFFFFFFFULL, 0x00000001FFFFFFFFULL, 0x00000003FFFFFFFFULL, 0x00000007FFFFFFFFULL,
        0x0000000FFFFFFFFFULL, 0x0000001FFFFFFFFFULL, 0x0000003FFFFFFFFFULL, 0x0000007FFFFFFFFFULL,
        0x000000FFFFFFFFFFULL, 0x000001FFFFFFFFFFULL, 0x000003FFFFFFFFFFULL, 0x000007FFFFFFFFFFULL,
        0x00000FFFFFFFFFFFULL, 0x00001FFFFFFFFFFFULL, 0x00003FFFFFFFFFFFULL, 0x00007FFFFFFFFFFFULL,
        0x0000FFFFFFFFFFFFULL, 0x0001FFFFFFFFFFFFULL, 0x0003FFFFFFFFFFFFULL, 0x0007FFFFFFFFFFFFULL,
        0x000FFFFFFFFFFFFFULL, 0x001FFFFFFFFFFFFFULL, 0x003FFFFFFFFFFFFFULL, 0x007FFFFFFFFFFFFFULL,
        0x00FFFFFFFFFFFFFFULL, 0x01FFFFFFFFFFFFFFULL, 0x03FFFFFFFFFFFFFFULL, 0x07FFFFFFFFFFFFFFULL,
        0x0FFFFFFFFFFFFFFFULL, 0x1FFFFFFFFFFFFFFFULL, 0x3FFFFFFFFFFFFFFFULL, 0x7FFFFFFFFFFFFFFFULL,
        0xFFFFFFFFFFFFFFFFULL};

    /** Lookup table for least significant set bit in a byte. */
    static constexpr uint8_t lt_lo[256] = {
        0x00, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00,
        0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00,
        0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00,
        0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x06, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00,
        0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00,
        0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00,
        0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00,
        0x01, 0x00, 0x07, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00,
        0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00,
        0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00,
        0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x06, 0x00, 0x01, 0x00, 0x02, 0x00,
        0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00,
        0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00,
        0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00,
        0x02, 0x00, 0x01, 0x00};

    /** Lookup table for select on bytes. */
    /*! Entry at idx = 256*j + i equals the position of the
     * (j+1)-th set bit in byte i. Positions lie in the range \f$[0..7]\f$.
     */
    static constexpr uint8_t lt_sel[256 * 8] = {
        0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2,
        0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0,
        1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
        0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
        2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3,
        0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
        1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,

        0, 0, 0, 1, 0, 2, 2, 1, 0, 3, 3, 1, 3, 2, 2, 1, 0, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 0, 5, 5, 1, 5,
        2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 0, 6, 6, 1, 6, 2, 2, 1, 6, 3,
        3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2,
        1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 0, 7, 7, 1, 7, 2, 2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1,
        4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4,
        3, 3, 1, 3, 2, 2, 1, 7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
        2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,

        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 3, 3, 2, 0, 0, 0, 4, 0, 4, 4, 2, 0, 4, 4, 3, 4, 3, 3, 2, 0, 0, 0, 5, 0,
        5, 5, 2, 0, 5, 5, 3, 5, 3, 3, 2, 0, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 0, 0, 0, 6, 0, 6, 6, 2, 0, 6,
        6, 3, 6, 3, 3, 2, 0, 6, 6, 4, 6, 4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 0, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3,
        2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 0, 0, 0, 7, 0, 7, 7, 2, 0, 7, 7, 3, 7, 3, 3, 2, 0, 7, 7, 4,
        7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 0, 7, 7, 5, 7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5,
        4, 4, 3, 4, 3, 3, 2, 0, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6, 4, 4, 2, 6, 4, 4, 3, 4, 3,
        3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 4, 4, 3, 0, 0, 0, 0, 0,
        0, 0, 5, 0, 0, 0, 5, 0, 5, 5, 3, 0, 0, 0, 5, 0, 5, 5, 4, 0, 5, 5, 4, 5, 4, 4, 3, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0,
        0, 6, 0, 6, 6, 3, 0, 0, 0, 6, 0, 6, 6, 4, 0, 6, 6, 4, 6, 4, 4, 3, 0, 0, 0, 6, 0, 6, 6, 5, 0, 6, 6, 5, 6, 5, 5,
        3, 0, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 7, 7, 3, 0, 0, 0, 7,
        0, 7, 7, 4, 0, 7, 7, 4, 7, 4, 4, 3, 0, 0, 0, 7, 0, 7, 7, 5, 0, 7, 7, 5, 7, 5, 5, 3, 0, 7, 7, 5, 7, 5, 5, 4, 7,
        5, 5, 4, 5, 4, 4, 3, 0, 0, 0, 7, 0, 7, 7, 6, 0, 7, 7, 6, 7, 6, 6, 3, 0, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4, 6, 4,
        4, 3, 0, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 5, 0, 5, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 6, 0, 6, 6, 4, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 6, 0, 6, 6,
        5, 0, 0, 0, 6, 0, 6, 6, 5, 0, 6, 6, 5, 6, 5, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0,
        0, 0, 0, 7, 0, 0, 0, 7, 0, 7, 7, 4, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 7, 7, 5, 0, 0, 0, 7, 0, 7, 7, 5, 0,
        7, 7, 5, 7, 5, 5, 4, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 7, 7, 6, 0, 0, 0, 7, 0, 7, 7, 6, 0, 7, 7, 6, 7, 6,
        6, 4, 0, 0, 0, 7, 0, 7, 7, 6, 0, 7, 7, 6, 7, 6, 6, 5, 0, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        6, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 6, 0, 6, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 7, 0,
        0, 0, 7, 0, 7, 7, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 7,
        7, 6, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 7, 7, 6, 0, 0, 0, 7, 0, 7, 7, 6, 0, 7, 7, 6, 7, 6, 6, 5,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 7, 0, 7, 7, 6,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7};

    /** Use to help to decide if a prefix sum stored in a byte overflows. */
    static constexpr uint64_t ps_overflow[65] = {
        0x8080808080808080ULL, 0x7f7f7f7f7f7f7f7fULL, 0x7e7e7e7e7e7e7e7eULL, 0x7d7d7d7d7d7d7d7dULL,
        0x7c7c7c7c7c7c7c7cULL, 0x7b7b7b7b7b7b7b7bULL, 0x7a7a7a7a7a7a7a7aULL, 0x7979797979797979ULL,
        0x7878787878787878ULL, 0x7777777777777777ULL, 0x7676767676767676ULL, 0x7575757575757575ULL,
        0x7474747474747474ULL, 0x7373737373737373ULL, 0x7272727272727272ULL, 0x7171717171717171ULL,
        0x7070707070707070ULL, 0x6f6f6f6f6f6f6f6fULL, 0x6e6e6e6e6e6e6e6eULL, 0x6d6d6d6d6d6d6d6dULL,
        0x6c6c6c6c6c6c6c6cULL, 0x6b6b6b6b6b6b6b6bULL, 0x6a6a6a6a6a6a6a6aULL, 0x6969696969696969ULL,
        0x6868686868686868ULL, 0x6767676767676767ULL, 0x6666666666666666ULL, 0x6565656565656565ULL,
        0x6464646464646464ULL, 0x6363636363636363ULL, 0x6262626262626262ULL, 0x6161616161616161ULL,
        0x6060606060606060ULL, 0x5f5f5f5f5f5f5f5fULL, 0x5e5e5e5e5e5e5e5eULL, 0x5d5d5d5d5d5d5d5dULL,
        0x5c5c5c5c5c5c5c5cULL, 0x5b5b5b5b5b5b5b5bULL, 0x5a5a5a5a5a5a5a5aULL, 0x5959595959595959ULL,
        0x5858585858585858ULL, 0x5757575757575757ULL, 0x5656565656565656ULL, 0x5555555555555555ULL,
        0x5454545454545454ULL, 0x5353535353535353ULL, 0x5252525252525252ULL, 0x5151515151515151ULL,
        0x5050505050505050ULL, 0x4f4f4f4f4f4f4f4fULL, 0x4e4e4e4e4e4e4e4eULL, 0x4d4d4d4d4d4d4d4dULL,
        0x4c4c4c4c4c4c4c4cULL, 0x4b4b4b4b4b4b4b4bULL, 0x4a4a4a4a4a4a4a4aULL, 0x4949494949494949ULL,
        0x4848484848484848ULL, 0x4747474747474747ULL, 0x4646464646464646ULL, 0x4545454545454545ULL,
        0x4444444444444444ULL, 0x4343434343434343ULL, 0x4242424242424242ULL, 0x4141414141414141ULL,
        0x4040404040404040ULL};

    /**
     *   @brief Counts the number of set bits in x.
     *
     *   @param  x 64-bit word
     *   @return Number of set bits.
     */
    static constexpr uint64_t cnt( uint64_t x );

    /**
     *   @brief Position of the most significant set bit the 64-bit word x
     *
     *   @param x 64-bit word
     *   @return The position (in 0..63) of the most significant set bit
     *   in `x` or 0 if x equals 0.
     *
     *   @sa sel, lo
     */
    static constexpr uint32_t hi( uint64_t x );

    /**
     *   @brief Calculates the position of the rightmost 1-bit in the 64bit
     *          integer x if it exists
     *
     *   @param x 64 bit integer.
     *
     *   @return The position (in 0..63) of the rightmost 1-bit in the 64bit
     *           integer x if x>0 and 0 if x equals 0.
     *
     *   @sa sel, hi
     */
    static constexpr uint32_t lo( uint64_t x );

    /**
     *   @brief Calculate the position of the i-th rightmost 1 bit in the 64bit integer x
     *
     *   @param x 64bit integer.
     *   @param i Argument i must be in the range \f$[1..cnt(x)]\f$.
     *
     *   @pre Argument i must be in the range \f$[1..cnt(x)]\f$.
     *   @sa hi, lo
     */
    static constexpr uint32_t sel( uint64_t x, uint32_t i );
    static constexpr uint32_t _sel( uint64_t x, uint32_t i );
  };

  // ============= inline - implementations ================

  // using built-in method or
  // 64-bit version of 32-bit proposal of
  // http://www-graphics.stanford.edu/~seander/bithacks.html
  template< typename T >
  constexpr uint32_t
  bits_impl< T >::hi( uint64_t x )
  {
#ifdef __SSE4_2__
    if ( x == 0 ) return 0;
    return 63 - __builtin_clzll( x );
#else
    uint64_t t{}, tt{};          // temporaries
    if ( ( tt = x >> 32 ) ) {    // hi >= 32
      if ( ( t = tt >> 16 ) ) {  // hi >= 48
        return ( tt = t >> 8 ) ? 56 + lt_hi[ tt ] : 48 + lt_hi[ t ];
      }
      else {  // hi < 48
        return ( t = tt >> 8 ) ? 40 + lt_hi[ t ] : 32 + lt_hi[ tt ];
      }
    }
    else {                      // hi < 32
      if ( ( t = x >> 16 ) ) {  // hi >= 16
        return ( tt = t >> 8 ) ? 24 + lt_hi[ tt ] : 16 + lt_hi[ t ];
      }
      else {  // hi < 16
        return ( tt = x >> 8 ) ? 8 + lt_hi[ tt ] : lt_hi[ x ];
      }
    }
#endif
  }

  // details see: http://citeseer.ist.psu.edu/leiserson98using.html
  // or page 10, Knuth TAOCP Vol 4 F1A
  template< typename T >
  constexpr uint32_t
  bits_impl< T >::lo(uint64_t x)
  {
#ifdef __SSE4_2__
    if ( x == 0 ) return 0;
    return __builtin_ctzll( x );
#else
    if ( x & 1 ) return 0;
    if ( x & 3 ) return 1;
    if ( x & 7 ) return 2;
    if ( x & 0x7FF ) {  // in average every second random number x can be
                        // answered this way
      return lt_lo[ ( x & 0x7FF ) >> 3 ] + 3;
    }
    // x&-x equals x with only the lsb set
    return lt_deBruijn_to_idx[ ( ( x & -x ) * deBruijn64 ) >> 58 ];
#endif
  }

  // see page 11, Knuth TAOCP Vol 4 F1A
  template< typename T >
  constexpr uint64_t
  bits_impl< T >::cnt( uint64_t x )
  {
#ifdef __SSE4_2__
    return __builtin_popcountll( x );
#else
#  ifdef POPCOUNT_TL
    return lt_cnt[ x & 0xFFULL ] + lt_cnt[ ( x >> 8 ) & 0xFFULL ]
           + lt_cnt[ ( x >> 16 ) & 0xFFULL ] + lt_cnt[ ( x >> 24 ) & 0xFFULL ]
           + lt_cnt[ ( x >> 32 ) & 0xFFULL ] + lt_cnt[ ( x >> 40 ) & 0xFFULL ]
           + lt_cnt[ ( x >> 48 ) & 0xFFULL ] + lt_cnt[ ( x >> 56 ) & 0xFFULL ];
#  else
    x = x - ( ( x >> 1 ) & 0x5555555555555555ull );
    x = ( x & 0x3333333333333333ull ) + ( ( x >> 2 ) & 0x3333333333333333ull );
    x = ( x + ( x >> 4 ) ) & 0x0f0f0f0f0f0f0f0full;
    return ( 0x0101010101010101ull * x >> 56 );
#  endif
#endif
  }

  template< typename T >
  constexpr uint32_t
  bits_impl< T >::sel( uint64_t x, uint32_t i )
  {
#if defined( __BMI__ ) && defined( __BMI2__ )
    // taken from folly
    return _tzcnt_u64( _pdep_u64( 1ULL << ( i - 1 ), x ) );
#endif
#ifdef __SSE4_2__
    uint64_t s = x, b{};
    s = s - ( ( s >> 1 ) & 0x5555555555555555ULL );
    s = ( s & 0x3333333333333333ULL ) + ( ( s >> 2 ) & 0x3333333333333333ULL );
    s = ( s + ( s >> 4 ) ) & 0x0F0F0F0F0F0F0F0FULL;
    s = 0x0101010101010101ULL * s;
    // now s contains 8 bytes s[7],...,s[0]; s[j] contains the cumulative sum
    // of (j+1)*8 least significant bits of s
    b = ( s + ps_overflow[ i ] ) & 0x8080808080808080ULL;
    // ps_overflow contains a bit mask x consisting of 8 bytes
    // x[7],...,x[0] and x[j] is set to 128-j
    // => a byte b[j] in b is >= 128 if cum sum >= j

    // __builtin_ctzll returns the number of trailing zeros, if b!=0
    int byte_nr = __builtin_ctzll( b ) >> 3;  // byte nr in [0..7]
    s <<= 8;
    i -= ( s >> ( byte_nr << 3 ) ) & 0xFFULL;
    return ( byte_nr << 3 )
           + lt_sel[ ( ( i - 1 ) << 8 )
                     + ( ( x >> ( byte_nr << 3 ) ) & 0xFFULL ) ];
#endif
    return _sel( x, i );
  }

  template< typename T >
  constexpr uint32_t
  bits_impl< T >::_sel( uint64_t x, uint32_t i )
  {
    uint64_t s = x, b{};  // s = sum
    s = s - ( ( s >> 1 ) & 0x5555555555555555ULL );
    s = ( s & 0x3333333333333333ULL ) + ( ( s >> 2 ) & 0x3333333333333333ULL );
    s = ( s + ( s >> 4 ) ) & 0x0F0F0F0F0F0F0F0FULL;
    s = 0x0101010101010101ULL * s;
    b = ( s + ps_overflow[ i ] );  //&0x8080808080808080ULL;// add something to
                                   //the partial sums to cause overflow
    i = ( i - 1 ) << 8;
    if ( b & 0x0000000080000000ULL )    // byte <=3
      if ( b & 0x0000000000008000ULL )  // byte <= 1
        if ( b & 0x0000000000000080ULL )
          return lt_sel[ ( x & 0xFFULL ) + i ];
        else
          return 8
                 + lt_sel[ ( ( ( x >> 8 ) & 0xFFULL ) + i
                             - ( ( s & 0xFFULL ) << 8 ) )
                           & 0x7FFULL ];  // byte 1;
      else                                // byte >1
        if ( b & 0x0000000000800000ULL )  // byte <=2
          return 16
                 + lt_sel[ ( ( ( x >> 16 ) & 0xFFULL ) + i
                             - ( s & 0xFF00ULL ) )
                           & 0x7FFULL ];  // byte 2;
        else
          return 24
                 + lt_sel[ ( ( ( x >> 24 ) & 0xFFULL ) + i
                             - ( ( s >> 8 ) & 0xFF00ULL ) )
                           & 0x7FFULL ];  // byte 3;
    else                                  //  byte > 3
      if ( b & 0x0000800000000000ULL )    // byte <=5
        if ( b & 0x0000008000000000ULL )  // byte <=4
          return 32
                 + lt_sel[ ( ( ( x >> 32 ) & 0xFFULL ) + i
                             - ( ( s >> 16 ) & 0xFF00ULL ) )
                           & 0x7FFULL ];  // byte 4;
        else
          return 40
                 + lt_sel[ ( ( ( x >> 40 ) & 0xFFULL ) + i
                             - ( ( s >> 24 ) & 0xFF00ULL ) )
                           & 0x7FFULL ];  // byte 5;
      else                                // byte >5
        if ( b & 0x0080000000000000ULL )  // byte<=6
          return 48
                 + lt_sel[ ( ( ( x >> 48 ) & 0xFFULL ) + i
                             - ( ( s >> 32 ) & 0xFF00ULL ) )
                           & 0x7FFULL ];  // byte 6;
        else
          return 56
                 + lt_sel[ ( ( ( x >> 56 ) & 0xFFULL ) + i
                             - ( ( s >> 40 ) & 0xFF00ULL ) )
                           & 0x7FFULL ];  // byte 7;
    return 0;
  }

  template< typename T >
  constexpr uint8_t bits_impl< T >::lt_cnt[ 256 ];
  template< typename T >
  constexpr uint32_t bits_impl< T >::lt_deBruijn_to_idx[ 64 ];
  template< typename T >
  constexpr uint32_t bits_impl< T >::lt_hi[ 256 ];
  template< typename T >
  constexpr uint64_t bits_impl< T >::lo_set[ 65 ];
  template< typename T >
  constexpr uint64_t bits_impl< T >::ps_overflow[ 65 ];
  template< typename T >
  constexpr uint8_t bits_impl< T >::lt_sel[ 256 * 8 ];
  template< typename T >
  constexpr uint8_t bits_impl< T >::lt_lo[ 256 ];

  using bits = bits_impl<>;

} // end namespace diverg

#endif  /* === #ifndef DIVERG_BITS_HPP__ === */
