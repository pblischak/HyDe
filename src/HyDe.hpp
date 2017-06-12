#ifndef HYDE_HPP
#define HYDE_HPP

/* Main class for running the HyDe test for hybridization detection.
 * Called from within main() with method run(). Parses command line options
 * and populates class member variables with input values. */

/* Include these here so that they available in main.cpp as well. */
#include <vector>

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
  void _printOut(const std::string& p1,
                 const std::string& hyb,
                 const std::string& p2,
                 const double& z, const double& p,
                 const double cp[16][16],
                 const double& avgObs,
                 const double& nObs,
                 std::ofstream& out);

  /* Functions for calculating stuff.*/
  inline int _convert(char str);                                                 /* Convert DNA bases from characters to ints.      */
  double _resolveAmbiguity(const int& out, const int& p1,                        /* Resolves ambiguity codes.                       */
                           const int& hyb, const int& p2,
                           double cp[16][16]);
  double _calcGH(const double cp[16][16], const double& nObs,                    /* Calculates GH test statistic.                   */
                 const double& avgObs, const int& p1,
                 const int& hyb, const int& p2);
  double _calcPvalue(const double& myZ);                                         /* Calculate p-value for GH test statistic.        */
  double _calcPvalueTwo(const double& myZ);
  void   _bootstrap();                                                           /* Bootstrap individuals within taxa.              */
  void   _resampleTaxonMap();                                                    /* Resample individuals within taxon map.          */
  inline long int _nCk(const long int& n, long int& k);                          /* Calculates binomial coefficient: n choose k.    */
  double _bonferroniCorrect();                                                   /* Bonferroni correction based on number of taxa.  */
  double _getCountMatrix(const int& p1,                                          /* Populates count matrix from given triplet.      */
                         const int& hyb,
                         const int& p2,
                         double cp[16][16]);

  /* private member variables */
  std::string _infile = "none", _mapfile = "none", _prefix = "hyde";
  double _pValue = 0.05;
  int _nInd = -999, _nSites = -999, _nTaxa = -999, _threads = 1, _bootReps = -999;
  std::vector<std::string> _indNames;
};

/* Convert DNA bases to ints as they are read in. */
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
    default : {std::cerr << "\n** ERROR: Unrecognized base character in DNA matrix: \"" << str << "\". **\n" << std::endl;
               exit(EXIT_FAILURE);} break;
  }

  return _baseCode;
}

/* Calculate binomial coefficient. */
inline long int HyDe::_nCk(const long int& n, long int& k){
  long int res = 1;
  if(k > n - k)
    k = n - k;

  for(long int i = 0; i < k; i++){
    res *= (n - i);
    res /= (i + 1.0);
  }
  return res;
}

#endif //HYDE_HPP
