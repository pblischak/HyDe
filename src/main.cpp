#include <iostream>
#include <string>
#include <cstring>

#include "HyDe.hpp"

/* Version info. */
const std::string Version = "0.1.0-alpha";
const std::string Date    = "February 2017";

void usage(); /* Print usage/helo information. */

int main(int argc, char* argv[]){
  /* Initial parsing of command line options to look for -h / --help,
   * -v / --version, or an incorrect number of arguments. */
  if(argc < 2){
    std::cerr << "\n** ERROR: Incorrect number of arguments. **" << std::endl;
    usage();
    exit(EXIT_FAILURE);
  }

  if(strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0){
    usage();
    exit(EXIT_SUCCESS);
  }

  if(strcmp(argv[1],"-v") == 0 || strcmp(argv[1],"--version") == 0){
    std::cerr << "\nThis is HyDe version " << Version << " (" << Date << ").\n" << std::endl;
    exit(EXIT_SUCCESS);
  }
  /* Pass argc and argv to constructor for HyDe class. It does a more intense
   * level of command line option parsing and initializes all variables. */
   //std::cout << "Please press ENTER to continue: ";
   //std::string in = "";
   //in = std::cin.get();
   HyDe hydeObj(argc, argv);
   hydeObj.run();

  return EXIT_SUCCESS;
}

void usage(){
  std::cerr << "\nUsage: hyde -i <infile> -m <taxon-map> -o <outgroup> -n <num-individuals> -t <num-taxa> -s <num-sites>  [additional options]" << std::endl
            << "\nInformation options:" << std::endl
            << "  -h [--help]          Prints this help message" << std::endl
            << "  -v [--version]       Prints version information" << std::endl
            << "\nProgram options:" << std::endl
            << "  -i [--infile]        Name of the data input file" << std::endl
            << "  -m [--map]           Map of individuals to taxa" << std::endl
            << "  -o [--outgroup]      Name of the outgroup (only one accepted)" << std::endl
            << "  -n [--num-ind]       Number of individuals in data matrix" << std::endl
            << "  -t [--num-taxa]      Number of taxa (species, OTUs)" << std::endl
            << "  -s [--num-sites]     Number of sites in the data matrix" << std::endl
            << "\nAdditional options:" << std::endl
            << "  -p [--pvalue]        P-value cutoff for test of significance (default=0.05)" << std::endl
            << "  -b [--bootstrap]     Number of bootstrap replicates to perform (default=0)" << std::endl
            << "  --prefix             Append a prefix to the beginning of outfile and logfile\n" << std::endl;
}
