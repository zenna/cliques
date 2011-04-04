// File: main_random_graph.cpp
// -- Generates a random graph using the benchmark approach 
//--------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks" (2008) V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
// and on the article "Dynamics and Modular Structure in Networks" (2008) R. Lambiotte, J.-C. Delvenne, M. Barahona
//
//--------------------------------------------------------------
// Copyright (c) 2010 Y. William Yu
//
// This program must not be distributed without agreement of the above author.
//--------------------------------------------------------------
// Author   : Y. William Yu
// Email    : yun.yu09@imperial.ac.uk
// Location : London, United Kingdom
// Time     : May 2010
//--------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <string>
#include <dlfcn.h>
#include <stdio.h>
#include <unistd.h>
#include <map>
#include <list>
#include <vector>


// Only used for development. Can be commented out later
#include <boost/date_time/posix_time/posix_time.hpp>

#include "graph_binary.h"

// Current version of the software
std::string version = "0.1";
namespace po = boost::program_options;

int
main(int argc, char **argv){
  // Following line is just for computing running times. Can be removed for actual code
  boost::posix_time::ptime t_start(boost::posix_time::microsec_clock::local_time());
  
  try {
    // Begin configuration parsing section
    std::string config_file;
    // Declare the command-line options
    po::options_description generic("Generic options");
    generic.add_options()
      ("help,h", "produce help message")
      ("version,v","print version string")
  //    ("config,c",po::value<std::string>(&config_file)->default_value("community.cfg"), "name of a file of a configuration.")
    ;

    // Declare the command-line and config file options
    po::options_description config("Configuration");
    config.add_options()
      ("weighted,w",po::value<bool>()->default_value(false)->zero_tokens(), "read the graph as a weighted one (weights are set to 1 otherwise)")
      ("random.n1", po::value<int>(), "Number of nodes in subgroup")
      ("random.n2", po::value<int>(), "Number of subgroups in group")
      ("random.n3", po::value<int>(), "Number of groups")
      ("random.k1", po::value<int>(), "Average subgroup node degree")
      ("random.k2", po::value<int>(), "Average non-subgroup, within group node degree")
      ("random.k3", po::value<int>(), "Average non-group node degree")
    ;
     
    // Declare hidden options, which are allowed on both command-line and config file, but don't show up on the help menu
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("graph-file", po::value<std::string>(),"output_file: write the generated random graph to this file (as binary)") 
    ;

    // Group together options_description(s)
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    // Pass all positional options to graph-file
    po::positional_options_description p;
    p.add("graph-file", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);
/*
    std::ifstream ifs(config_file.c_str());
    if (!ifs) {
      std::cout << "can not open config file: " << config_file << '\n';
      return 0;
    }
    else {
      store(parse_config_file(ifs, config_file_options),vm);
      po::notify(vm);
    }
*/
    if (vm.count("help")||(argc==1)){ 
      std::cout << "Usage: " << argv[0] << " output_file [options]" << '\n';
      std::cout << visible << '\n';
      return 1;
    }
    if (vm.count("version")) {
      std::cout << "Version " << version << '\n';
      return 1;
    }

    // End parsing section

    // Initialize random number generator with current time
    std::srand(time(NULL));

    Graph g(vm["random.n1"].as<int>(),vm["random.k1"].as<int>(),vm["random.n2"].as<int>(),vm["random.k2"].as<int>(),vm["random.n3"].as<int>(),vm["random.k3"].as<int>());

    g.display_binary(vm["graph-file"].as<std::string>().c_str());
    g.display();
  }
  catch(std::exception& e)
  {
    std::cout << e.what() << '\n';
    return 1;
  }



  //boost::posix_time::ptime t_end(boost::posix_time::microsec_clock::local_time());
  //std::cout << "Time elapsed: " << t_end-t_start << '\n';

  //std::cout << "test";

  return 0;
}

