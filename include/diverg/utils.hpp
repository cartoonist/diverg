/**
 *    @file  utils.hpp
 *   @brief  Utility and helper functions.
 *
 *  This header file contains general utility and helper functions.
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

#ifndef  DIVERG_UTILS_HPP__
#define  DIVERG_UTILS_HPP__

#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cctype>
#include <set>

#include <sdsl/enc_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <gum/iterators.hpp>


#define DIVERG_DEFAULT_TMPDIR "/tmp"
#define DIVERG_TMPFILE_TEMPLATE "/diverg-XXXXXX"

#define BINARY_NAME "diverg"

#ifdef __ASSERT_VOID_CAST
#define DIVERG_ASSERT_VOID_CAST __ASSERT_VOID_CAST
#else
#define DIVERG_ASSERT_VOID_CAST static_cast<void>
#endif  // ifdef __ASSERT_VOID_CAST

#define DIVERG_ASSERT(expr)							\
  ((expr)								\
   ? DIVERG_ASSERT_VOID_CAST (0)						\
   : diverg::assert_fail (#expr, BINARY_NAME, __FILE__, __LINE__, __PRETTY_FUNCTION__))

#define DIVERG_MACRO_MIN(x, y) ((x) < (y) ? (x) : (y))
#define DIVERG_MACRO_MAX(x, y) ((x) < (y) ? (y) : (x))


namespace diverg {
  /**
   *  @brief  Assert fail function
   *
   *  @param  expr The assert expression.
   *  @param  outfile The name of output file.
   *  @param  file The source file that assertion failed.
   *  @param  line The line number where assertion failed.
   *  @param  func The function signiture where assertion failed.
   *
   *  Print the message and exit with error code 134.
   */
    inline void
  assert_fail( std::string const& expr, std::string const& outfile,
      std::string const& file, int line, std::string const& func )
  {
    std::cout << outfile << ": " << file << ":" << line << ": " << func
      << ": Assertion `" << expr << "' failed." << "\n"
      << "Aborted." << std::endl;
    std::exit( 134 );
  }


  /**
   *  @brief  Check whether a string ends with another string.
   *
   *  @param  str The first string.
   *  @param  suf The second string.
   *  @return `true` if `suf` is a suffix of `str`; otherwise `false`.
   *
   *  It checks the first string whether the second one is one of its suffixes or not.
   *
   *  Specialized for `std::string` and `std::string`.
   *
   *  NOTE: `suf` is passed by VALUE!
   */
    inline bool
  ends_with( const std::string& str, std::string suf )
  {
    if ( suf.length() <= str.length() &&
        std::equal( suf.rbegin(), suf.rend(), str.rbegin() ) ) {
      return true;
    }
    return false;
  }  /* -----  end of function ends_with  ----- */


  /**
   *  @brief  Check whether a string starts with another string.
   *
   *  @param  str The first string.
   *  @param  pre The second string.
   *  @return `true` if `pre` is a prefix of `str`; otherwise `false`.
   *
   *  It checks the first string whether the second one is one of its prefixes or not.
   *
   *  Specialized for `std::string` and `std::string`.
   *
   *  NOTE: `pre` is passed by VALUE!
   */
    inline bool
  starts_with( const std::string& str, std::string pre )
  {
    if ( pre.length() <= str.length() &&
        std::equal( pre.begin(), pre.end(), str.begin() ) ) {
      return true;
    }
    return false;
  }  /* -----  end of function starts_with  ----- */


    inline std::string
  get_env( const std::string& var )
  {
    const char* val = ::getenv( var.c_str() );
    if ( val == 0 ) {
      return "";
    }
    else {
      return val;
    }
  }

    inline std::string
  get_tmpdir_env( )
  {
    return get_env( "TMPDIR" );
  }

    inline std::string
  get_tmpdir( )
  {
    std::string tmpdir = get_tmpdir_env();
    if ( tmpdir.size() == 0 ) tmpdir = DIVERG_DEFAULT_TMPDIR;
    return tmpdir;
  }

    inline std::string
  get_tmpfile( char const* directory="" )
  {
    assert( std::strlen( directory ) == 0 || starts_with( directory, "/" ) );
    std::string tfpath = get_tmpdir() + directory + DIVERG_TMPFILE_TEMPLATE;
    char* tmpl = new char [ tfpath.size() + 1 ];
    std::strcpy( tmpl, tfpath.c_str() );
    int fd = mkstemp( tmpl );
    tfpath = tmpl;

    ::close( fd );
    delete[] tmpl;
    return tfpath;
  }


  typedef uint64_t DefaultContainerSize;

  /**
   *  @brief  Simple object serialization implementation.
   *
   *  @param[out]  out The output stream.
   *  @param[in]  obj The object to be written into the output stream.
   *
   *  It simply writes the object bitwise into the stream.
   */
  template< typename TObject >
      inline void
    serialize( std::ostream& out, const TObject& obj )
    {
      out.write( reinterpret_cast< const char* >( &obj ), sizeof( TObject ) );
    }

  /**
   *  @brief  Container serialization implementation.
   *
   *  @param[out]  out The output stream.
   *  @param[in]  size The size of the container.
   *  @param[in]  begin The begin input iterator.
   *  @param[in]  end The end input iterator.
   *
   *  First, the size of the container is written represented by `TSize` type followed
   *  by all items between two input iterators are written.
   */
  template< typename TIter, typename TSize >
    inline void
  _serialize( std::ostream& out, TSize size, TIter begin, TIter end )
  {
    serialize( out, size );
    for ( ; begin != end; ++begin ) serialize( out, *begin );
  }  /* -----  end of template function _serialize  ----- */

  /**
   *  @brief  Serialize an ordered container to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  begin The begin input iterator.
   *  @param[in]  end The end input iterator.
   *
   *  It gets two iterators of a container and serialize it to the given output stream
   *  by writing each item to the stream following the size of the container. The
   *  `TSize` represents the size type.
   */
  template< typename TIter, typename TSize = DefaultContainerSize >
    inline void
  serialize( std::ostream& out, TIter begin, TIter end )
  {
    TSize size = end - begin;
    _serialize( out, size, begin, end );
  }  /* -----  end of template function serialize  ----- */

  /**
   *  @brief  Serialize an unordered container to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  container The container.
   *  @param[in]  begin The begin input iterator.
   *  @param[in]  end The end input iterator.
   *
   *  It gets two iterators of a container and serialize it to the given output stream
   *  by writing each item to the stream following the size of the container. The
   *  `TSize` represents the size type. Since the container is unordered, it needs
   *  container itself to infer the size.
   */
  template< typename TContainer, typename TIter, typename TSize = DefaultContainerSize >
    inline void
  serialize( std::ostream& out, const TContainer& container, TIter begin, TIter end )
  {
    TSize size = container.size();
    _serialize( out, size, begin, end );
  }  /* -----  end of template function serialize  ----- */

  /**
   *  @brief  Serialize an `deque` to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  v The `deque`.
   *  @param[in]  begin The begin input iterator [unused].
   *  @param[in]  end The end input iterator [unused].
   *
   *  A wraper to `serialize` member function of `deque`.
   */
  template< typename T >
    inline void
  serialize( std::ostream& out, const std::deque< T >& v )
  {
    serialize( out, v.begin(), v.end() );
  }  /* -----  end of template function serialize  ----- */

  /**
   *  @brief  Serialize an `vector` to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  v The `vector`.
   *  @param[in]  begin The begin input iterator [unused].
   *  @param[in]  end The end input iterator [unused].
   *
   *  A wraper to `serialize` member function of `vector`.
   */
  template< typename T >
    inline void
  serialize( std::ostream& out, const std::vector< T >& v )
  {
    serialize( out, v.begin(), v.end() );
  }  /* -----  end of template function serialize  ----- */

  /**
   *  @brief  Serialize an `enc_vector` to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  ev The `enc_vector`.
   *  @param[in]  begin The begin input iterator [unused].
   *  @param[in]  end The end input iterator [unused].
   *
   *  A wraper to `serialize` member function of `enc_vector`.
   */
  template< typename TCoder, uint32_t TDens = 128, uint8_t TWidth=0 >
    inline void
  serialize( std::ostream& out, const sdsl::enc_vector< TCoder, TDens, TWidth >& ev )
  {
    ev.serialize( out );
  }  /* -----  end of template function serialize  ----- */

  /**
   *  @brief  Serialize an `int_vector_buffer` to an output stream.
   *
   *  @param[out]  out Output stream.
   *  @param[in]  ivb The `int_vector_buffer`.
   */
  template< uint8_t TWidth >
  inline void
  serialize( std::ostream& out, const sdsl::int_vector_buffer< TWidth >& ivb )
  {
    throw std::runtime_error( "`sdsl::int_vector_buffer` cannot be serialised" );
  }  /* -----  end of template function serialize  ----- */


  /**
   *  @brief  Deserialize a simple object from an input stream.
   *
   *  @param[in,out]  in The input stream.
   *  @param[out]  obj The deserialized value will be written to this variable.
   *
   *  It simply reads the object from an input stream.
   */
  template< typename TObject >
      inline void
    deserialize( std::istream& in, TObject& obj )
    {
      in.read( reinterpret_cast< char* >( &obj ), sizeof( TObject ) );
    }


  template< typename TObject, typename ...TArgs >
      inline void
    open( TObject& obj, std::istream& in, TArgs&&... args )
    {
      obj.load( in, std::forward< TArgs >( args )... );
    }

  template< typename TObject, typename ...TArgs >
      inline void
    open( TObject& obj, const std::string& file_name, TArgs&&... args )
    {
      std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
      if( !ifs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      open( obj, ifs, std::forward< TArgs >( args )... );
    }  /* -----  end of template function load  ----- */


  template< typename TObject, typename ...TArgs >
      inline void
    save( TObject& obj, std::ostream& out, TArgs&&... args )
    {
      obj.serialize( out, std::forward< TArgs >( args )... );
    }

  template< typename TObject, typename ...TArgs >
      inline void
    save( TObject& obj, const std::string& file_name, TArgs&&... args )
    {
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      if( !ofs ) {
        throw std::runtime_error( "cannot open file '" + file_name + "'" );
      }
      save( obj, ofs, std::forward< TArgs >( args )... );
    }


  /**
   *  @brief  Clear the given container.
   *
   *  @param  c The container to be cleaned.
   *
   *  A wrapper function to provide an interface to `clean` member function of
   *  `std::vector`.
   */
  template< typename T >
      inline void
    clear( std::vector< T >& c )
    {
      c.clear();
    }

  /**
   *  @brief  Clear the given container.
   *
   *  @param  c The container to be cleaned.
   *
   *  A wrapper function to provide an interface to `clean` member function of
   *  `std::vector`.
   */
  template< typename T >
      inline void
    clear( std::deque< T >& c )
    {
      c.clear();
    }

  /**
   *  @brief  Clear the given container.
   *
   *  @param  c The container to be cleaned.
   *
   *  A wrapper function to provide an interface to `clean` member function of
   *  `sdsl::enc_vector`.
   */
  template< typename TCoder, uint32_t TDens = 128, uint8_t TWidth=0 >
      inline void
    clear( sdsl::enc_vector< TCoder, TDens, TWidth >& ev )
    {
      sdsl::util::clear( ev );
    }

  /**
   *  @brief  Clear the given container.
   *
   *  @param  ivb The container to be cleaned.
   *
   */
  template< uint8_t TWidth >
  inline void
  clear( sdsl::int_vector_buffer< TWidth >& ivb )
  {
    ivb.reset();
  }


  /**
   *  @brief  Reserve the required memory for the vector of size `size`.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  It calls reserve member function of the container.
   */
  template< typename TObject, typename TSize >
      inline void
    reserve( std::vector< TObject >& container, const TSize size )
    {
      container.reserve( size );
    }

  /**
   *  @brief  Reserve the required memory for the string of size `size`.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  It calls reserve member function of the container.
   */
  template< typename TSize >
      inline void
    reserve( std::string& container, const TSize size )
    {
      container.reserve( size );
    }

  /**
   *  @overload The `std::deque` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  Do nothing.
   */
  template< typename TObject, typename TSize >
      inline void
    reserve( std::deque< TObject >&, const TSize )
    {
      /* NOOP */
    }

  /**
   *  @overload The `std::set` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  Do nothing.
   */
  template< typename TObject, typename TSize >
      inline void
    reserve( std::set< TObject >&, const TSize )
    {
      /* NOOP */
    }

  /**
   *  @overload The `sdsl::enc_vector` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  Do nothing.
   */
  template< typename TCoder, typename TSize, uint32_t TDens = 128, uint8_t TWidth=0 >
      inline void
    reserve( sdsl::enc_vector< TCoder, TDens, TWidth >&, const TSize )
    {
      /* NOOP */
    }

  /**
   *  @overload The `sdsl::int_vector_buffer` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   *  Do nothing.
   */
  template< uint8_t TWidth, typename TSize >
  inline void
  reserve( sdsl::int_vector_buffer< TWidth >&, const TSize )
  {
    /* NOOP */
  }

  /**
   *  @overload The `sdsl::int_vector` cannot be reserved.
   *
   *  @param  container The container.
   *  @param  size Size to reserve.
   *
   */
  template< uint8_t TWidth, typename TSize >
  inline void
  reserve( sdsl::int_vector< TWidth >&, const TSize )
  {
    /* NOOP */
  }


  /**
   *  @brief The shrink to fit interface function.
   *
   *  It acts a wrapper for `shrink_to_fit` member function.
   */
  template< typename TContainer >
  inline void
  shrink_to_fit( TContainer& container )
  {
    container.shrink_to_fit();
  }

  /**
   *  @overload The shrink to fit interface function.
   */
  template< typename TCoder, uint32_t TDens = 128, uint8_t TWidth=0 >
  inline void
  shrink_to_fit( sdsl::enc_vector< TCoder, TDens, TWidth >& )
  {
    /* NOOP */
  }

  /**
   *  @overload The shrink to fit interface function.
   */
  template< uint8_t TWidth >
  inline void
  shrink_to_fit( sdsl::int_vector_buffer< TWidth >& )
  {
    /* NOOP */
  }


  /**
   *  @brief  Resize the container to the given `size`.
   *
   *  @param  container The container.
   *  @param  size Size to resize.
   *
   *  It calls resize member function of the container.
   */
  template< typename TContainer, typename TSize >
  inline void
  resize( TContainer& container, const TSize size )
  {
    container.resize( size );
  }

  /**
   *  @overload The `sdsl::int_vector_buffer` cannot be resized.
   *
   *  @param  container The container.
   *  @param  size Size to resize.
   *
   *  Do nothing.
   */
  template< typename TCoder, typename TSize, uint32_t TDens = 128, uint8_t TWidth=0 >
  inline void
  resize( sdsl::enc_vector< TCoder, TDens, TWidth >&, const TSize )
  {
    /* NOOP */
  }

  /**
   *  @overload The `sdsl::int_vector_buffer` cannot be resized.
   *
   *  @param  container The container.
   *  @param  size Size to resize.
   *
   *  Do nothing.
   */
  template< uint8_t TWidth, typename TSize >
  inline void
  resize( sdsl::int_vector_buffer< TWidth >& ivb, const TSize size )
  {
    ivb[ size - 1 ] = 0;
  }


  /**
   *  @brief  Deserialize a container from an input stream.
   *
   *  @param[in,out]  in The input stream.
   *  @param[out]  container The output container.
   *  @param[out]  itr The output iterator.
   *
   *  It gets an output container, and an output iterator. Then reads from input stream
   *  and populate the container from serialized data.
   */
  template< typename TContainer, typename TOutIter, typename TSize = DefaultContainerSize >
    inline void
  deserialize( std::istream& in, TContainer& container, TOutIter itr )
  {
    TSize size;
    deserialize( in, size );
    reserve( container, size );
    for ( unsigned int i = 0; i < size; ++i ) {
      typename TContainer::value_type item;
      deserialize( in, item );
      *itr++ = std::move( item );
    }
  }  /* -----  end of template function deserialize  ----- */

  template< typename T >
    inline void
  deserialize( std::istream& in, std::vector< T >& container )
  {
    deserialize( in, container, std::back_inserter( container ) );
  }  /* -----  end of template function deserialize  ----- */

  template< typename T >
    inline void
  deserialize( std::istream& in, std::deque< T >& container )
  {
    deserialize( in, container, std::back_inserter( container ) );
  }  /* -----  end of template function deserialize  ----- */

  template< typename TCoder, uint32_t TDens = 128, uint8_t TWidth=0 >
    inline void
  deserialize( std::istream& in, sdsl::enc_vector< TCoder, TDens, TWidth >& ev )
  {
    ev.load( in );
  }  /* -----  end of template function deserialize  ----- */

  template< uint8_t TWidth >
    inline void
  deserialize( std::istream& in, sdsl::int_vector_buffer< TWidth >& ivb )
  {
    throw std::runtime_error( "`sdsl::int_vector_buffer` cannot be deserialised" );
  }  /* -----  end of template function serialize  ----- */
}  /* --- end of namespace diverg --- */

#endif  /* --- #ifndef DIVERG_UTILS_HPP__ --- */
