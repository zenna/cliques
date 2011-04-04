// File: main_community.cpp
// -- community detection 
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

#include "abstract_community_finder.h"

// constant string codes to make coloring the output easier; probably a bit of overkill,  but all in good fun.
const std::string black_text= "\033[22;30m";
const std::string red_text = "\033[22;31m";
const std::string green_text = "\033[22;32m";
const std::string brown_text = "\033[22;33m";
const std::string blue_text= "\033[22;34m";
const std::string magenta_text = "\033[22;35m";
const std::string cyan_text = "\033[22;36m";
const std::string white_text = "\033[22;37m";
const std::string default_text = "\033[0m";

// Current version of the software
std::string version = "0.0.2";
namespace po = boost::program_options;

// Below couple of lines are for reading dynamic libraries
// size of buffer for reading in directory entries
static unsigned int BUF_SIZE = 1024;
// our global factory for making AbstractCommunityFinders
std::map<std::string, maker_t *, std::less<std::string> > factory;

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
      ("config,c",po::value<std::string>(&config_file)->default_value("community.cfg"), "name of a file of a configuration.")
    ;

    // Declare the command-line and config file options
    po::options_description config("Configuration");
    config.add_options()
      ("weighted,w",po::value<bool>()->default_value(false)->zero_tokens(), "read the graph as a weighted one (weights are set to 1 otherwise)")
      ("precision,q", po::value<double>(), "set epsilon -- a given pass stops when the quality is increased by less than epsilon")
      ("timet,t", po::value<double>(), "Resolution parameter. Optimisation of stability in order to uncover smaller communities. Modularity: t=1. t must be smaller than 1. Alternately, \"Markov time\"")
      ("display,l", po::value<int>()->default_value(-2), "Displays the graph of level arg rather than the hierarchical structure.")
      ("verbose", po::value<bool>()->default_value(false)->zero_tokens(), "Display config options. Disable for production runs.")
    ;
     
    // Declare hidden options, which are allowed on both command-line and config file, but don't show up on the help menu
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("graph-file", po::value<std::string>(),"input_file: read the graph to partition from this file") 
      ("number-runs", po::value<unsigned int>(),"number of times to run the partition optimization")
      ("save-intermediates",po::value<bool>()->default_value(false)->zero_tokens(), "choose whether or not to save intermediate runs of the partition optimization") // Unused!!
      ("math.quality-function", po::value<std::string>(),"choose the quality function to use")
      ("math.partition-method", po::value<std::string>(),"choose the partition method to use")
      ("math.community-type", po::value<std::string>(),"choose community library to load")
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

    std::ifstream ifs(config_file.c_str());
    if (!ifs) {
      std::cout << "can not open config file: " << config_file << '\n';
      return 0;
    }
    else {
      store(parse_config_file(ifs, config_file_options),vm);
      po::notify(vm);
    }



    if (vm.count("help")||(argc==1)){ 
      std::cout << "Usage: " << argv[0] << " input_file [options]" << '\n';
      std::cout << visible << '\n';
      return 1;
    }
    if (vm.count("version")) {
      std::cout << "Version " << version << '\n';
      return 1;
    }
    
    if (vm["verbose"].as<bool>()){
      std::cout << "\033[22;35mCommunity finder v" << version << default_text << '\n';
      if (vm["weighted"].as<bool>()){
        std::cout << "Reading " << green_text << "weighted" << default_text <<" graph ..." << '\n';
      }
      else std::cout << "Reading " << green_text << "unweighted" << default_text << " graph ..." << '\n';

      if (vm.count("precision")){std::cout << green_text << "Precision" << default_text << " = " << vm["precision"].as<double>() << '\n';}
        else {std::cout << "Please specify a " << red_text << " precision " << default_text << " level" << '\n'; return 1;}
      if (vm.count("timet")){std::cout << "Markov" << green_text << " timet " << default_text << " = " << vm["timet"].as<double>() << '\n';}
        else {std::cout << "Please specify " << red_text << " timet"<< default_text << "!" << '\n'; return 1;}
      {std::cout << green_text << "Display" << default_text << " level = "<< vm["display"].as<int>() << '\n';}
      if (vm.count("graph-file")) std::cout << "Input file is: " << vm["graph-file"].as<std::string>() << '\n';
    }
    // End parsing section

    // Begin dynamic library loading section. *.so in the lib directory will be loaded. All of the *.so files have makers associated that will automatically register with the factory.
    FILE *dl; // handle to read directory
    char *command_str = "ls lib/*.so"; // command string to get dynamic lib names
    char in_buf[BUF_SIZE]; // input buffer for lib names
    std::list<void *> dl_list; // list to hold handles for dynamic libs
    std::list<void *>::iterator itr;
    std::vector<std::string> AbstractCommunityFinder_names; // vector of AbstractCommunityFinder types to build list
    std::list<AbstractCommunityFinder *> AbstractCommunityFinder_list; // list of AbstractCommunityFinder objects we create
    std::list<AbstractCommunityFinder *>::iterator sitr;
    std::map<std::string, maker_t *, std::less<std::string> >::iterator fitr;
    // get the names of all the dynamic libs (.so files) in the current dir
    dl = popen(command_str, "r");
    if(!dl){
      perror("popen");
      exit(-1);
    }
    void *dlib;
    char name[1024];
    while (fgets(in_buf, BUF_SIZE, dl)){
      // trim off the whitespace
      char *ws = strpbrk(in_buf, " \t\n");
      if(ws) *ws = '\0';
      // append ./ to the front of the lib name
      sprintf(name, "./%s", in_buf);
      dlib = dlopen(name, RTLD_NOW);
      if(dlib == NULL){
        std::cerr << dlerror() << '\n';
        exit(-1);
      }
      // add the handle to our list
      dl_list.insert(dl_list.end(), dlib);
    }
    int i = 0;
    // create an array of the AbstractCommunityFinder names
    for (fitr=factory.begin(); fitr!=factory.end(); fitr++){
      AbstractCommunityFinder_names.insert(AbstractCommunityFinder_names.end(),fitr->first);
      i++;
    }

    // Initialize random number generator with current time
    std::srand(time(NULL));

    Graph g(vm["graph-file"].as<std::string>().c_str() , vm["weighted"].as<bool>());

    int display_level = vm["display"].as<int>();
    
    int finder_choice = -1;
    if (vm["verbose"].as<bool>()) std::cout << "Possible choices for type of community:" << '\n';
    for (unsigned int i=0; i<AbstractCommunityFinder_names.size(); i++) {
      if (vm["verbose"].as<bool>()) std::cout << '\t' << "(" << i << ") " << AbstractCommunityFinder_names[i] << '\n';
      if (AbstractCommunityFinder_names[i] == vm["math.community-type"].as<std::string>()) finder_choice = i;
    }
    if (finder_choice!=-1) {
      if (vm["verbose"].as<bool>()) std::cout << "Running using community type: " << green_text << AbstractCommunityFinder_names[finder_choice] << default_text << '\n';
    }
    else {
      std::cerr << "Community choice "<<red_text << vm["math.community-type"].as<std::string>() << default_text<< " did not match any available. Please make certain you have all the relevant libraries properly installed and compiled." << '\n';
      return(1);
    }

    // Insert string buffers
    std::ostringstream stdout_current;
    std::ostringstream stderr_current;
    // Make sure that we save at least one partition
    double quality_current= -1000000000;
    // Create the Abstract Community Finder object
    AbstractCommunityFinder_list.insert(AbstractCommunityFinder_list.end(),factory[AbstractCommunityFinder_names[finder_choice]](g, -1, vm["precision"].as<double>()));
    // Run the find_partition algorithm the proscribed number of times and keep the partition of the highest quality.
    sitr=AbstractCommunityFinder_list.begin();
    for (unsigned int i=0; i<vm["number-runs"].as<unsigned int>(); i++){
      AbstractCommunityFinder::result findings=(*sitr)->find_partition(vm["timet"].as<double>(),display_level);
      if (findings.quality > quality_current) {
        quality_current = findings.quality;
        stdout_current.str("");
        stderr_current.str("");
        stdout_current << findings.partition;
        stderr_current << findings.partition_summary;
      }
    }

    std::cout << stdout_current.str();
    std::cerr << stderr_current.str();

    // Delete all objects allocated, and release the dynamic libraries
    for(sitr=AbstractCommunityFinder_list.begin();sitr!=AbstractCommunityFinder_list.end();sitr++){
      delete *sitr;
    }
    for(itr=dl_list.begin(); itr!=dl_list.end() ; itr++){
      dlclose(*itr);
    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what() << '\n';
    return 1;
  }



  boost::posix_time::ptime t_end(boost::posix_time::microsec_clock::local_time());
  //std::cout << "Time elapsed: " << t_end-t_start << '\n';

  //std::cout << "test";

  return 0;
}

