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


template< typename TOrdinal = int32_t, typename TSize = uint64_t >
struct Options {
  /* === TYPE MEMBERS === */
  using ordinal_type = TOrdinal;
  using size_type = TSize;
  /* === DATA MEMBERS === */
  std::string graph_path;
  std::string format;
  TOrdinal n;
  TSize nnz;
  unsigned int seed;
  unsigned int dlo;
  unsigned int dup;
  std::string seg_name;
  std::string reg_name;
  std::string partition;
  std::string grid;
  unsigned int l1size;
  bool run_kokkos;
  bool run_rspgemm;
  bool host;
  std::string out_path;
  bool compare;
  int verbose;
  /* === LIFECYCLE === */
  Options()
      : graph_path( "" ), format( "" ), n( 0 ), nnz( 0 ), dlo( 0 ), dup( 0 ),
        seg_name( "" ), reg_name( "" ), partition( "team-sequential" ),
        grid( "auto" ), l1size( 8192 ), run_kokkos( true ),
        run_rspgemm( true ), host( false ), out_path( "" ), compare( false ),
        verbose( 0 )
  { }
};

template< typename TOrdinal = Options<>::ordinal_type,
          typename TSize = Options<>::size_type >
Options< TOrdinal, TSize >
parse_arguments( int argc, char* argv[] )
{
  Options< TOrdinal, TSize > opts;

  // Parse command line arguments
  for ( int i = 1; i < argc; i++ ) {
    if ( ( strcmp( argv[ i ], "-g" ) == 0 ) || ( strcmp( argv[ i ], "--graph" ) == 0 ) ) {
      opts.graph_path = argv[ ++i ];
      std::cout << "Parameter: graph <- " << opts.graph_path << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-f" ) == 0 ) || ( strcmp( argv[ i ], "--format" ) == 0 ) ) {
      opts.format = argv[ ++i ];
      std::cout << "Parameter: format <- " << opts.format << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-n" ) == 0 ) || ( strcmp( argv[ i ], "--order" ) == 0 ) ) {
      opts.n = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "Parameter: n <- " << opts.n << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-z" ) == 0 ) || ( strcmp( argv[ i ], "--nnz" ) == 0 ) ) {
      opts.nnz = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "Parameter: nnz <- " << opts.nnz << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-S" ) == 0 ) || ( strcmp( argv[ i ], "--seed" ) == 0 ) ) {
      opts.seed = atoi( argv[ ++i ] );
      std::cout << "Parameter: seed <- " << opts.seed << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-d" ) == 0 ) || ( strcmp( argv[ i ], "--dist-lower" ) == 0 ) ) {
      opts.dlo = atoi( argv[ ++i ] );
      std::cout << "Parameter: d <- " << opts.dlo << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-D" ) == 0 ) || ( strcmp( argv[ i ], "--dist-upper" ) == 0 ) ) {
      opts.dup = atoi( argv[ ++i ] );
      std::cout << "Parameter: D <- " << opts.dup << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-K" ) == 0 ) || ( strcmp( argv[ i ], "--disable-kokkos" ) == 0 ) ) {
      opts.run_kokkos = false;
      std::cout << "Parameter: disable-kokkos <- "
                << ( opts.run_kokkos ? "false" : "true" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-R" ) == 0 ) || ( strcmp( argv[ i ], "--disable-rspgemm" ) == 0 ) ) {
      opts.run_rspgemm = false;
      std::cout << "Parameter: disable-rspgemm <- "
                << ( opts.run_rspgemm ? "false" : "true" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-s" ) == 0 ) || ( strcmp( argv[ i ], "--segment-name" ) == 0 ) ) {
      opts.seg_name = argv[ ++i ];
      std::cout << "Parameter: segment-name <- " << opts.seg_name << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-r" ) == 0 ) || ( strcmp( argv[ i ], "--region-name" ) == 0 ) ) {
      opts.reg_name = argv[ ++i ];
      std::cout << "Parameter: region-name <- " << opts.reg_name << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-p" ) == 0 ) || ( strcmp( argv[ i ], "--partition" ) == 0 ) ) {
      opts.partition = argv[ ++i ];
      std::cout << "Parameter: partition <- " << opts.partition << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-G" ) == 0 ) || ( strcmp( argv[ i ], "--grid" ) == 0 ) ) {
      opts.grid = argv[ ++i ];
      std::cout << "Parameter: grid <- " << opts.grid << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-l" ) == 0 ) || ( strcmp( argv[ i ], "--l1-size" ) == 0 ) ) {
      opts.l1size = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "Parameter: l1size <- " << opts.l1size << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-H" ) == 0 ) || ( strcmp( argv[ i ], "--host" ) == 0 ) ) {
      opts.host = true;
      std::cout << "Parameter: host <- "
                << ( opts.host ? "true" : "false" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-o" ) == 0 ) || ( strcmp( argv[ i ], "--output" ) == 0 ) ) {
      opts.out_path = argv[ ++i ];
      std::cout << "Parameter: output <- " << opts.graph_path << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-c" ) == 0 ) || ( strcmp( argv[ i ], "--compare" ) == 0 ) ) {
      opts.compare = true;
      std::cout << "Parameter: compare <- "
                << ( opts.compare ? "true" : "false" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-v" ) == 0 ) || ( strcmp( argv[ i ], "--verbose" ) == 0 ) ) {
      opts.verbose++;
    }
    else if ( ( strcmp( argv[ i ], "-h" ) == 0 ) || ( strcmp( argv[ i ], "--help" ) == 0 ) ) {
      std::cout
          << "Options:\n"
          << "  * Adjacency matrix as input (overrides -n and -z)\n"
          << "    --graph (-g) <path>:       graph file path\n"
          << "    --format (-f) <str>:       specifies graph file format when it is not\n"
          << "                               automatically detected; e.g. when GFA has no header\n"
          << "                               (formats: 'gfa', 'gfa1', 'gfa2')\n"
          << "    --region-name (-r) <str>:  specifies target region by path id (embedded in graph)\n"
          << "                               (overrides -s)\n"
          << "    --segment-name (-s) <str>: specifies target region by the name of its first node\n"
          << "                               or segment (default: first node/segment in the graph)\n"
          << "  * Random matrix as input\n"
          << "    --order (-n) <int>:        exponent num, determines number of rows 2^num\n"
          << "                               (default: 2^12 = 4096 or 2^((z+1)/2) if -z given)\n"
          << "    --nnz (-z) <int>:          exponent num, determines total matrix size 2^num\n"
          << "                               (default: 2^22 = 4096*1024 or 2^n*1024 if -n provided)\n"
          << "    --rng-seed (-S) <int>:     Seed for random number generation (default: random)\n"
          << "  --dist-lower (-d) <int>:     lower bound of outer distance, i.e. fragment size\n"
          << "                               (default: 0)\n"
          << "  --dist-upper (-D) <int>:     upper bound of outer distance, i.e. fragment size\n"
          << "                               (required)\n"
          << "  --partition (-p) <str>:      Partitioning schemes: 'team-sequential' (default),\n"
          << "                               'thread-sequential', 'thread-parallel', and\n"
          << "                               'thread-range-paralell'\n"
          << "  --grid (-G) <str>            Execution grid options: 'auto' (auto-tuning; default),\n"
          << "                               'suggested' (based on row density of input)\n"
          << "  --grid (-G) <int,int[,int]>  '<team size>,<vector size>[,<team work size>]'\n"
          << "                               in which default team work size is "
          << RCRS_BENCHMARK_DEFAULT_TEAM_WORK_SIZE_STR << "\n"
          << "                               (any non-numeric, non-whitespace delimiter works)\n"
          << "                               examples: '-G 32,4,16' or '-G 16x16'\n"
          << "  --l1-size (-l) <int>:        exponent num, determines L1 size of accumulator 2^num;\n"
          << "                               with constraints: 10<=num<=15 (default: 2^13 = 1KB)\n"
          << "  --host (-H):                 run on host (default: false)\n"
          << "  --disable-kokkos (-K):       disable benchmarking with Kokkos SpGEMM (PairG)\n"
          << "                               (default: false)\n"
          << "  --disable-rspgemm (-R):      disable benchmarking with rSpGEMM (DiVerG)\n"
          << "                               (default: false)\n"
          << "  --output (-o) <path>:        specifies output path where DiVerG index will be saved\n"
          << "  --compare (-c):              compare resulting indices build by DiVerG and PairG\n"
          << "                               given that -K is not set (default: false)\n"
          << "  --verbose (-v):              increase verbose level by one (default: 0)\n"
          << "  --help (-h):                 print this message"
          << std::endl;
      exit( EXIT_SUCCESS );
    }
    else {
      std::cerr << "[ERROR] Unknown argument '" << argv[ i ] << "'" << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  if ( opts.verbose ) std::cout << "Verbose level: " << opts.verbose << std::endl;

  return opts;
}

template< typename TOrdinal, typename TSize >
void
check_options( Options< TOrdinal, TSize >& opts )
{
  if ( opts.graph_path.empty() ) {
    if ( opts.nnz == 0 ) {
      if ( opts.n == 0 )  // If both opts.n and opts.nnz are not provided:
        opts.nnz = pow( 2, 22 );
      else  // If opts.n is provided but opts.nnz is not:
        opts.nnz = opts.n * 1024;
      std::cout << "[DEFAULT] nnz = " << opts.nnz << std::endl;
    }

    // If opts.n is not provided, set order to 2^((log2(opts.nnz)+1)/2).
    if ( opts.n == 0 ) {
      int power = 1;
      auto temp = opts.nnz;
      while ( ( temp /= 2 ) ) ++power;
      opts.n = pow( 2, power / 2 );
      std::cout << "[DEFAULT] n = " << opts.n << std::endl;
    }

    if ( opts.n < 0 /* || opts.nnz < 0 */ ) {
      std::cerr << "[ERROR] nnz must be greater than zero" << std::endl;
      exit( EXIT_FAILURE );
    }

    if ( opts.nnz / opts.n > static_cast< TSize >( opts.n ) ) {
      std::cerr << "[ERROR] Non-zero values cannot be fit (nnz > n*n)"
                << std::endl;
      exit( EXIT_FAILURE );
    }
  }
  else {
    if ( !gum::util::readable( opts.graph_path ) ) {
      std::cerr << "[ERROR] Input file is missing or not readable"
                << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  if ( opts.dup == 0 ) {
    std::cerr << "[ERROR] Too few arguments (required argument: '-D')"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  if ( opts.dup < opts.dlo ) {
    std::cerr << "[WARN] Lower bound of distance constraints is larger than "
                 "upper bound"
              << std::endl;
    std::swap( opts.dlo, opts.dup );
    std::cerr << "[WARN] Swapped distance constraints bounds: (" << opts.dlo
              << ", " << opts.dup << ")" << std::endl;
  }

  if ( !opts.run_kokkos && !opts.run_rspgemm ) {
    std::cerr << "[ERROR] Nothing to do!" << std::endl;
    exit( EXIT_FAILURE );
  }
}

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
