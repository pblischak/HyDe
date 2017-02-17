#ifndef HYDE_HPP
#define HYDE_HPP

/* Main class for running the HyDe test for hybridization detection.
 * Called from within main() with method run(). Parses command line options
 * and populates class member variables with input values. */

/* Include these here so that they available in main.cpp as well. */
 #include <vector>
 #include <unordered_map>

class HyDe {
public:
  HyDe(int c, char* v[]);
  ~HyDe(){};
  void run();

private:
  /* private member functions */
  /* Functions for reading and writing information. */
  void _parseCommandLine(int ac, char* av[]);
  void _checkCommandLineInput();
  void _parseSpeciesMap();
  void _readInfile();
  void _writeOoutfile();
  void _writeLogFile();

  /* Functions for calculating .*/
  inline int _convert(char str);          /* Convert DNA bases form characters to ints.      */
  void _resolveAmbiguity(int code);        /* Randomly resolves ambiguity codes.              */
  double _getGH();                         /* Calculates GH test statistic.                   */
  double _getPvalue();                     /* Calculate p-value with bonferroni correction.   */
  double _bonferroniCorrect();             /* Bonferroni correctiono based on number of taxa. */

  /* private member variables */
  std::string _infile = "none", _mapfile = "none", _logfile = "log.txt",
              _outfile = "out.txt", _prefix = "hyde", _outgroup = "none";
  double _pValue = 0.05;
  int _nInd = -999, _nSites = -999, _nTaxa = -999;
  std::vector<std::string> _indNames, _speciesNames;
  std::unordered_map<std::string, std::vector<int> > _speciesMap;
  std::vector<int>  _dnaMatrix, _uniqueSitePatterns, _uniqueSitePatternCounts;
  bool _log = 0, _quiet = 0;

};

inline int HyDe::_convert(char str){
  int _baseCode = 999;
  switch (str) {
    case 'A': case 'a': _baseCode = 0;  break;
    case 'G': case 'g': _baseCode = 1;  break;
    case 'C': case 'c': _baseCode = 2;  break;
    case 'T': case 't': case 'U': case 'u': _baseCode = 3;  break;
    case 'M': case 'm': _baseCode = 5;  break;
    case 'R': case 'r': _baseCode = 6;  break;
    case 'W': case 'w': _baseCode = 7;  break;
    case 'S': case 's': _baseCode = 8;  break;
    case 'Y': case 'y': _baseCode = 9;  break;
    case 'K': case 'k': _baseCode = 10;  break;
    case 'B': case 'b': _baseCode = 11;  break;
    case 'D': case 'd': _baseCode = 12;  break;
    case 'H': case 'h': _baseCode = 13;  break;
    case 'V': case 'v': _baseCode = 14;  break;
    case 'N': case 'n': _baseCode = 15;  break;
    case '-': _baseCode = 4; break;
    case '?': _baseCode = 15; break;
    default : break;
  }

  return _baseCode;
}

#endif //HYDE_HPP
