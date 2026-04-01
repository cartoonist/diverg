/**
 *    @file  eti_macros.hpp
 *   @brief  Macros for Explicit Template Instantiation (ETI) support.
 *
 *  When DIVERG_ETI_ENABLED is defined (compiled library in use),
 *  DIVERG_ETI_EXTERN_CLASS suppresses implicit instantiation by declaring
 *  the class specialisation extern. The compiled library provides the
 *  one and only instantiation.
 *
 *  When DIVERG_ETI_ENABLED is not defined (header-only mode), both macros
 *  expand to nothing and implicit instantiation proceeds as normal.
 *
 *  Usage in a header (after the class definition, outside the class):
 *    DIVERG_ETI_EXTERN_CLASS( MyClass<int, double> )
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Tue Mar 31, 2026  00:00
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2026, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_ETI_MACROS_HPP__
#define DIVERG_ETI_MACROS_HPP__

#ifdef DIVERG_ETI_ENABLED
#  define DIVERG_ETI_EXTERN_CLASS(...)  extern template class __VA_ARGS__;
#else
#  define DIVERG_ETI_EXTERN_CLASS(...)
#endif

#endif  /* --- #ifndef DIVERG_ETI_MACROS_HPP__ --- */
