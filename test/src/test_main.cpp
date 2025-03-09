/**
 *    @file  test_main.cpp
 *   @brief  Unit test framework (Catch2) implementation file.
 *
 *  This file is the custom main entry point for the unit tests. It initialises
 *  and finalies the Kokkos runtime. It also registers a custom event listener
 *  to set the internal random generator seed from the command line argument.
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

#include <iostream>

#include<Kokkos_Core.hpp>

#include "catch2/catch_session.hpp"
#include "catch2/catch_get_random_seed.hpp"
#include "catch2/reporters/catch_reporter_event_listener.hpp"
#include "catch2/reporters/catch_reporter_registrars.hpp"

#include "test_base.hpp"


class MyEventListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void
    set_rnd_seed()
    {
      auto seed = Catch::getSeed();
      if ( seed != 0 ) {
        std::cout << "Setting random generator seed to " << seed << "..."
                  << std::endl;
        rnd::set_seed( seed );
      }
    }

    void
    testRunStarting( Catch::TestRunInfo const& ) override
    {
      this->set_rnd_seed();
    }
};

CATCH_REGISTER_LISTENER( MyEventListener )

int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );

  int result = Catch::Session().run( argc, argv );

  Kokkos::finalize();

  return result;
}
