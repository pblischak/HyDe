#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>

#include "HyDe.hpp"

/*************************************************/
/* Public member functions accessed from main(). */
/*************************************************/

/* Constructor. */
HyDe::HyDe(int c, char* v[]){
  _parseCommandLine(c, v);
  _checkCommandLineInput();
  _parseSpeciesMap();
  _dnaMatrix.resize(_nInd * _nSites);
  //_readInfile();
}

/* Runs HyDe analysis. */
void HyDe::run(){
  std::cerr << "We're running!" << std::endl;
}

/**************************************************************/
/* Private member functions accessed from within HyDe class . */
/**************************************************************/

/* Read in and parse command line flags. */
void HyDe::_parseCommandLine(int ac, char* av[]){
  int _invalidArgCount = 0;
  std::vector<std::string> _invalidArgs;
  for(int i = 1; i < ac; i++){
    if(strcmp(av[i], "-i") == 0 || strcmp(av[i], "--infile") == 0){
      _infile = av[i + 1];
    } else if(strcmp(av[i], "-n") == 0 || strcmp(av[i], "--num-ind") == 0){
      _nInd = atoi(av[i + 1]);
    } else if(strcmp(av[i], "-s") == 0 || strcmp(av[i], "--num-sites") == 0){
      _nSites = atoi(av[i + 1]);
    } else if(strcmp(av[i], "-t") == 0 || strcmp(av[i], "--num-taxa") == 0){
      _nTaxa = atoi(av[i + 1]);
    } else if(strcmp(av[i], "-m") == 0 || strcmp(av[i], "--map") == 0){
      _mapfile = av[i + 1];
    } else if(strcmp(av[i], "-o") == 0 || strcmp(av[i], "--outgroup") == 0){
      _outgroup = av[i + 1];
    } else if(strcmp(av[i], "-p") == 0 || strcmp(av[i], "--pvalue") == 0){
      _pValue = atof(av[i + 1]);
    } else if(strcmp(av[i], "--prefix") == 0){
      _prefix = av[i + 1];
    } else if(strcmp(av[i], "-q") == 0 || strcmp(av[i], "--quiet") == 0){
      _quiet = 1;
    } else if(strcmp(av[i], "--log") == 0){
      _log = 1;
    } else if(av[i][0] == '-'){ /* This checks if it is a flag that isn't valid. */
      _invalidArgCount++;
      _invalidArgs.push_back(av[i]);
    }
  }

  if(_invalidArgCount != 0){
    std::cerr << "\n** ERROR: Unrecognized command line flag(s). **\n" << std::endl;
    for(int i = 0; i < _invalidArgs.size(); i++){
      std::cerr << "  " << _invalidArgs[i] << std::endl;
    }
    std::cerr << "\nType 'hyde -h' for command line options.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
}

/* Check that valid aptions were passed to all parameters. */
void HyDe::_checkCommandLineInput(){
  int _errorCaught = 0;
  if(!_quiet)
    std::cerr << "\nChecking command line input...            ";

  if(_infile.compare("none") == 0){
    std::cerr << "\nMissing or invalid option for -i [--infile].\n";
    _errorCaught++;
  }
  if(_nInd < 0){
    std::cerr << "\nMissing or invalid option for -n [--num-ind].\n";
    _errorCaught++;
  }
  if(_nSites < 0){
    std::cerr << "\nMissing or invalid option for -s [--num-sites].\n";
    _errorCaught++;
  }
  if(_nTaxa < 0){
    std::cerr << "\nMissing or invalid option for -t [--num-taxa].\n";
    _errorCaught++;
  }
  if(_mapfile.compare("none") == 0){
    std::cerr << "\nMissing or invalid option for -m [--map].\n";
    _errorCaught++;
  }
  if(_outgroup.compare("none") == 0){
    std::cerr << "\nMissing or invalid option for -o [--outgroup].\n";
    _errorCaught++;
  }

  if(_errorCaught > 0){
    std::cerr << "** ERROR: " << _errorCaught << " command line options were improperly specified. **\n" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    std::cerr << "Good" << std::endl;
  }
}

/* Read in species map file. */
void HyDe::_parseSpeciesMap(){
  std::ifstream _mapStream(_mapfile);
  std::string _s1, _s2, _s2prev = "";
  int _indCount = 0, _taxaCount = 0;
  if(!_mapStream.is_open()){
    while(_mapStream >> _s1 >> _s2){
      _indCount++;
      _speciesMap.push_back(_s1);
      _speciesMap.push_back(_s2);
      if(_s2.compare(_s2prev) != 0){
        _taxaCount++;
        _s2prev = _s2;
      }
    }
  } else {
    std::cerr << "\n** ERROR: Cannot open species map file. **\n" << std::endl
              << "  File specified: " << _mapfile << "\n" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(_indCount != _nInd){
    std::cerr << "\n** ERROR: Number of individuals in species map does not match the number specified at the command line. **\n" << std::endl
              << "  Number from command line:   " << _nInd << std::endl
              << "  Number in species map file: " << _indCount << std::endl << std::endl;
  }

  if(_taxaCount != _nTaxa){
    std::cerr << "\n** ERROR: Number of taxa in species map does not match the number specified at the command line. **\n" << std::endl
              << "  Number from command line:   " << _nTaxa << std::endl
              << "  Number in species map file: " << _taxaCount << std::endl << std::endl;
  }
}

void HyDe::_readInfile(){
  std::ifstream _infileStream(_infile);
  std::string _str1, _str2;
  int _indIndex = 0, _errCount = 0;
  if(!_infileStream.is_open()){

  } else {
    while(_infileStream >> _str1 >> _str2){
      if(_str1.compare(_speciesMap[_indIndex * 2]) != 0){
        std::cerr << "\n** ERROR: Name in infile does not match name in map file. **\n" << std::endl
                  << "  Line " << _indIndex + 1 << ": " << _str1 << " vs." << _speciesMap[_indIndex * 2] << "\n" << std::endl;
        _errCount++;
      }

      for(int s = 0; s < _str2.length(); s++){
        _dnaMatrix[_indIndex * _nSites + s] = _convert(_str2[s]);
      }
      _indIndex++;
    }
    if(_errCount > 0){
      exit(EXIT_FAILURE);
    }
  }

}
