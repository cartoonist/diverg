/**
 *    @file  random.hpp
 *   @brief  Random number generation module.
 *
 *  This header file contains general utility and helper functions for random
 *  number generation.
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

#ifndef DIVERG_RANDOM_HPP__
#define DIVERG_RANDOM_HPP__

#include <random>
#include <cassert>
#include <limits>


namespace diverg {
  namespace random {
    /* Adapted from here: https://stackoverflow.com/q/440133/357257 */
    thread_local static std::random_device rd;
    thread_local static std::mt19937 gen( rd() );

    template< typename TFloat, typename TGenerator >
    inline TFloat
    random_real( TFloat low, TFloat high, TGenerator&& rgen )
    {
      assert( low < high );  // half-open range: [low, high) imposed by std::uniform_real_distribution<>
      std::uniform_real_distribution< TFloat > dis( low, high );
      return dis( rgen );
    }

    template< typename TFloat >
    inline TFloat
    random_real( TFloat low=0, TFloat high=1 )
    {
      return random_real( low, high, gen );
    }

    template< typename TInteger, typename TGenerator >
    inline TInteger
    random_integer( TInteger low,
                    TInteger high,
                    TGenerator&& rgen )
    {
      assert( low <= high );  // closed range: [low, high] imposed by std::uniform_int_distribution<>
      std::uniform_int_distribution< TInteger > dis( low, high );
      return dis( rgen );
    }

    template< typename TInteger >
    inline TInteger
    random_integer( TInteger low=std::numeric_limits< TInteger >::min(),
                    TInteger high=std::numeric_limits< TInteger >::max() )
    {
      return random_integer( low, high, gen );
    }

    template< typename TInteger, typename TGenerator >
    inline TInteger
    random_index( TInteger length, TGenerator&& rgen )
    {
      assert( 0 < length );
      return random_integer< TInteger >( 0, length - 1, rgen );
    }

    template< typename TInteger >
    inline TInteger
    random_index( TInteger length )
    {
      return random_index( length, gen );
    }

    template< typename TInteger, typename TGenerator >
    inline std::string
    random_string( TInteger length, TGenerator&& rgen,
                   const char* charset=nullptr, std::size_t charset_len=0 )
    {
      const char alphanum_charset[] = "0123456789"
                                      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                      "abcdefghijklmnopqrstuvwxyz";
      if ( charset == nullptr ) {
        charset = alphanum_charset;
        charset_len = sizeof( alphanum_charset ) - 1 /* null-terminated */;
      }

      auto randchar = [&rgen, &charset, &charset_len]() -> char
      {
        return charset[ random_index( charset_len, rgen ) ];
      };
      std::string str( length, 0 );
      std::generate_n( str.begin(), length, randchar );
      return str;
    }

    template< typename TInteger >
    inline std::string
    random_string( TInteger length,
                   const char* charset=nullptr, std::size_t charset_len=0 )
    {
      return random_string( length, gen, charset, charset_len );
    }
  }  /* --- end of namespace random --- */
}  /* --- end of namespace diverg --- */

#endif // DIVERG_RANDOM_HPP__
