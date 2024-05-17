/**
 *    @file  rcrs_benchmark.hpp
 *   @brief  Header file for benchmark Range CRS matrix operations
 *
 *  This header file provides the helper functions and utilities for
 *  benchmarking Range CRS; particularly for command-line argument parsing.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Mon May 13, 2024  11:51
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2024, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_RCRS_BENCHMARK_HPP__
#define DIVERG_RCRS_BENCHMARK_HPP__

#include <cmath>

#include <diverg/basic_types.hpp>
#include <diverg/dindex.hpp>

#define RCRS_BENCHMARK_DEFAULT_TEAM_WORK_SIZE 16
#define RCRS_BENCHMARK_DEFAULT_TEAM_WORK_SIZE_STR                             \
  std::to_string( RCRS_BENCHMARK_DEFAULT_TEAM_WORK_SIZE )


template< typename TGraph >
TGraph
load_graph( std::string const& graph_path, std::string const& format )
{
  TGraph graph;

  auto load_versioned_gfa = []( auto& graph, auto& graph_path, bool sorted,
                                auto version ) {
    using dynamic_type =
        typename std::remove_reference_t< decltype( graph ) >::dynamic_type;
    std::ifstream ifs( graph_path, std::ifstream::in | std::ifstream::binary );
    gfak::GFAKluge gg;
    gg.set_version( version );
    std::cout << "Parsing input graph (GFAKluge)..." << std::endl;
    gg.parse_gfa_file( ifs );
    dynamic_type dyn_graph;
    std::cout << "Constructing a Dynamic graph (gum)..." << std::endl;
    gum::util::extend( dyn_graph, gg, sorted );
    std::cout << "Converting to a Succinct graph (gum)..." << std::endl;
    graph = dyn_graph;
  };

  if ( format == "gfa" ) {
    std::cout << "Auto-detecting GFA version..." << std::endl;
    gum::util::load_gfa( graph, graph_path, true );
  }
  else if ( format == "gfa1" ) {
    std::cout << "Enforcing GFA version 1.0..." << std::endl;
    load_versioned_gfa( graph, graph_path, true, 1.0 );
  }
  else if ( format == "gfa2" ) {
    std::cout << "Enforcing GFA version 2.0..." << std::endl;
    load_versioned_gfa( graph, graph_path, true, 2.0 );
  }
  else if ( format == "" ) {
    std::cout << "Auto-detecting graph file format..." << std::endl;
    gum::util::load( graph, graph_path, true );
  }
  else
    throw std::runtime_error( "[ERROR] unknown file format '" + format + "'" );

  return graph;
}

template< typename TGraph >
auto
region_nodes_rank_range( TGraph const& graph, std::string const& reg_name,
                         std::string const& seg_name )
{
  using rank_type = typename TGraph::rank_type;

  std::pair< rank_type, rank_type > rank_range{ 0, 0 };

  if ( !reg_name.empty() ) {
    auto id = diverg::util::first_node_in_region( graph, reg_name );

    if ( id != 0 ) {
      rank_range.first = graph.id_to_rank( id );
      std::cout << "Found region '" << reg_name << "' with start node '"
                << graph.get_node_prop( rank_range.first ).name
                << "' with rank " << rank_range.first << std::endl;
    }
    else {
      std::cerr
          << "[WARN] There is no paths in graph corresponding to region '"
          << reg_name << "'" << std::endl;
    }
  }

  if ( rank_range.first == 0 && !seg_name.empty() ) {
    auto id = diverg::util::node_by_name( graph, seg_name );

    if ( id != 0 ) {
      rank_range.first = graph.id_to_rank( id );
      std::cout << "Found segment '" << seg_name << "' with rank "
                << rank_range.first << std::endl;
    }
    else {
      std::cerr << "[WARN] Segment '" << seg_name << "' not found"
                << std::endl;
    }
  }

  if ( rank_range.first == 0 ) rank_range.first = 1;

  rank_type nof_comps = 0;
  rank_type comp_i = 0;
  std::tie( rank_range.second, comp_i, nof_comps )
      = diverg::util::upper_node_rank_in_component( graph, rank_range.first );

  DIVERG_ASSERT( rank_range.first < rank_range.second );

  std::cout << "Found " << nof_comps << " components" << std::endl;
  std::cout << "Selecting component " << comp_i << " with "
            << rank_range.second - rank_range.first << " nodes ["
            << rank_range.first << ", " << rank_range.second << ")..."
            << std::endl;

  return rank_range;
}

#endif // DIVERG_RCRS_BENCHMARK_HPP__
