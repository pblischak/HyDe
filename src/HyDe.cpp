#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "HyDe.hpp"

/* Declaring "global" variables to be passed to threads.
 * They're not really global though because they are
 * only accessible to functions in 'HyDe.cpp'. In an earlier approach
 * all of these variables were members of the HyDe class, which required
 * copying them before passing to threads. I think that this prevents
 * needing to copy things. */
 double _counts1[16][16] = {0.0},
        _counts2[16][16] = {0.0},
        _counts3[16][16] = {0.0},
        _pVals[3] = {0.0},
        _zVals[3] = {0.0};
 int _numObs[3] = {0};
 std::vector<int> _dnaMatrix;
 std::unordered_map<std::string, std::vector<int> > _taxaMap;
 std::vector<std::string> _taxaNames;
 std::string _outgroup = "none";
 std::unordered_map<int, std::vector<int> > _baseLookup = { /* {A,G,C,T} == {0,1,2,3} */
  {0, {0}},
  {1, {1}},
  {2, {2}},
  {3, {3}},
  {4, {4}},           /* - == ignore */
  {5,  {0, 2}},       /* M == A or C */
  {6,  {0, 1}},       /* R == A or G */
  {7,  {0, 3}},       /* W == A ot T */
  {8,  {1, 2}},       /* S == G or C */
  {9,  {2, 3}},       /* Y == C or T */
  {10, {1, 3}},       /* K == G or T */
  {11, {1, 2, 3}},    /* B == C or G or T */
  {12, {0, 1, 3}},    /* D == A or G or T */
  {13, {0, 2, 3}},    /* H == A or C or T */
  {14, {0, 1, 2}},    /* V == A or G or C */
  {15, {0, 1, 2, 3}}  /* N == A or G or C or T */
 };

/*************************************************/
/* Public member functions accessed from main(). */
/*************************************************/

/* Constructor. */
HyDe::HyDe(int c, char* v[]){
  _parseCommandLine(c, v);
  _checkCommandLineInput();
  _parseSpeciesMap();
  _dnaMatrix.resize(_nInd * _nSites);
  _readInfile();
}

/* Runs HyDe analysis. */
void HyDe::run(){
  /* Open outfile and logfile (if specified). */
  std::string _outfile = _prefix + "-out.txt";
  std::ofstream _outStream;
  _outStream.open(_outfile, std::ios::out | std::ios::app);
  if(!_outStream.is_open()){
    std::cerr << "** ERROR: Could not open outfile: " << _outfile << ". **\n" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    _outStream << "P1\tHybrid\tP2\tZscore\tPvalue\tX1\tX2\tX3\tX4\tX5\tX6\tX7\tX8\tX9\tX10\tX11\tX12\tX13\tX14\tX15" << std::endl;
  }

  #pragma omp parallel for num_threads(_threads) schedule(dynamic) firstprivate(_counts1, _counts2, _counts3, _taxaMap, _taxaNames, _dnaMatrix, _baseLookup, _outgroup, _numObs, _zVals, _pVals)
  for(unsigned i = 0; i < _taxaNames.size() - 2; i++){
    for(unsigned j = i + 1; j < _taxaNames.size() - 1; j++){
      for(unsigned k = j + 1; k < _taxaNames.size(); k++){
        /* Re-initialize count matrices to 0. */
        for(int a = 0; a < 16; a++){
          for(int b = 0; b < 16; b++){
            _counts1[a][b] = 0.0;
            _counts2[a][b] = 0.0;
            _counts3[a][b] = 0.0;
          }
        }
        _numObs[0] = _getCountMatrix(_taxaNames[i], _taxaNames[j], _taxaNames[k], _counts1);
        _numObs[1] = _getCountMatrix(_taxaNames[j], _taxaNames[i], _taxaNames[k], _counts2);
        _numObs[2] = _getCountMatrix(_taxaNames[i], _taxaNames[k], _taxaNames[j], _counts3);
        //std::cerr << _numObs[0] << "\t" << _numObs[1] << "\t" << _numObs[2] << std::endl;

        /* Calculate the GH statistic. */
        _zVals[0] = _calcGH(_counts1, _numObs[0]);
        _zVals[1] = _calcGH(_counts2, _numObs[1]);
        _zVals[2] = _calcGH(_counts3, _numObs[2]);

        /* Get p-values. */
        _pVals[0] = _calcPvalue(_zVals[0]);
        _pVals[1] = _calcPvalue(_zVals[1]);
        _pVals[2] = _calcPvalue(_zVals[2]);

        #pragma omp critical
        {
          if(_zVals[0] != -99999.9)
            _printOut(_taxaNames[i], _taxaNames[j], _taxaNames[k], _zVals[0], _pVals[0], _counts1, _outStream);
          if(_zVals[1] != -99999.9)
            _printOut(_taxaNames[j], _taxaNames[i], _taxaNames[k], _zVals[1], _pVals[1], _counts2, _outStream);
          if(_zVals[2] != -99999.9)
            _printOut(_taxaNames[i], _taxaNames[k], _taxaNames[j], _zVals[2], _pVals[2], _counts3, _outStream);
        }
      }
    }
  }
}

/*************************************************************/
/* Private member functions accessed from within HyDe class. */
/*************************************************************/

/* Read in and parse command line flags. */
void HyDe::_parseCommandLine(int ac, char* av[]){
  int _invalidArgCount = 0;
  std::vector<char*> _invalidArgs;
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
    } else if(strcmp(av[i], "--threads") == 0){
      #ifdef _OPENMP
        _threads = atoi(av[i + 1]);
        int _maxThreads = omp_get_max_threads();
        if(_threads > _maxThreads){
          std::cerr << "\n** WARNING: Requesting more threads (" << _threads << ") than are available (" << _maxThreads << "). **" <<std::endl
                    << "   Defaulting to maximum number of threads available.\n" << std::endl;
          _threads = _maxThreads;
        }
      #else
        std::cerr << "\n** WARNING: OpenMP not enabled. --threads argument will be ignored. **\n" << std::endl;
      #endif
    } else if(strcmp(av[i], "--prefix") == 0){
      _prefix = av[i + 1];
    } else if(av[i][0] == '-'){ /* This checks if it is a flag that isn't valid. */
      _invalidArgCount++;
      _invalidArgs.push_back(av[i]);
    }
  }

  if(_invalidArgCount != 0){
    std::cerr << "\n** ERROR: Unrecognized command line flag(s). **\n" << std::endl;
    for(unsigned i = 0; i < _invalidArgs.size(); i++){
      std::cerr << "  " << _invalidArgs[i] << std::endl;
    }
    std::cerr << "\nType 'hyde -h' for command line options.\n" << std::endl;
    exit(EXIT_FAILURE);
  }
}

/* Check that valid options were passed to all parameters. */
void HyDe::_checkCommandLineInput(){
  int _errorCaught = 0;
  if(strcmp(_infile.c_str(), "none") == 0){
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
  if(strcmp(_mapfile.c_str(),"none") == 0){
    std::cerr << "\nMissing or invalid option for -m [--map].\n";
    _errorCaught++;
  }
  if(strcmp(_outgroup.c_str(), "none") == 0){
    std::cerr << "\nMissing or invalid option for -o [--outgroup].\n";
    _errorCaught++;
  }

  if(_errorCaught > 0){
    std::cerr << "** ERROR: " << _errorCaught << " command line option(s) were improperly specified. **\n" << std::endl;
    exit(EXIT_FAILURE);
  }
}

/* Read in taxon map file. */
void HyDe::_parseSpeciesMap(){
  std::ifstream _mapStream(_mapfile);
  std::string _str1, _str2, _str2prev = "";
  int _indCount = 0, _taxaCount = 0;
  bool _outgroupFound = 0;
  if(_mapStream.is_open()){
    while(_mapStream >> _str1 >> _str2){
      //std::cerr << _str1 << "\t" << _str2 << std::endl;
      _indNames.push_back(_str1);
      _taxaMap[_str2].push_back(_indCount);
      _indCount++;
      if(_str2.compare(_str2prev) != 0){
        _taxaCount++;
        if(_str2.compare(_outgroup) != 0){
          _taxaNames.push_back(_str2);
        } else {
          if(!_outgroupFound) _outgroupFound = 1;
        }
        _str2prev = _str2;
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
              exit(EXIT_FAILURE);
  }

  if(_taxaCount != _nTaxa){
    std::cerr << "\n** ERROR: Number of taxa in species map does not match the number specified at the command line. **\n" << std::endl
              << "  Number from command line:   " << _nTaxa << std::endl
              << "  Number in species map file: " << _taxaCount << std::endl << std::endl;
              exit(EXIT_FAILURE);
  }

  if(!_outgroupFound){
    std::cerr << "\n** ERROR: Outgroup specified at the command line (" << _outgroup << ") does not match any taxon names in the species map file.\n" << std::endl
              << "Taxon names in species map file:" << std::endl;
    for(unsigned i = 0; i < _taxaNames.size(); i++){
      std::cerr << "  " << _taxaNames[i] << std::endl;
    }
    exit(EXIT_FAILURE);
  }
}

/* Read in DNA matrix in sequential Phylip format w/o header info.        */
/* Converts DNA bases to integer codes using inlined _convert function.   */
void HyDe::_readInfile(){
  std::ifstream _infileStream(_infile);
  std::string _str1, _str2;
  int _indIndex = 0, _errCount = 0;
  if(!_infileStream.is_open()){
    std::cerr << "\n** ERROR: Cannot open specified infile: " << _infile << ". **\n" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    while(_infileStream >> _str1 >> _str2){
      if(_str1.compare(_indNames[_indIndex]) != 0){
        std::cerr << "\n** ERROR: Name in infile does not match name in map file. **\n" << std::endl
                  << "  Line " << _indIndex + 1 << ": " << _str1 << " vs." << _indNames[_indIndex] << "\n" << std::endl;
        _errCount++;
      } else if(_str2.length() != (unsigned) _nSites){
        std::cerr << "\n** ERROR: Length of input sequence not equal to specified number of sites (" << _nSites << ")\n" << std::endl
                  << "  Line " << _indIndex + 1 << " in " << _infile << " has " << _str2.length() << " sites.\n" << std::endl;
        _errCount++;
      }

      for(unsigned s = 0; s < _str2.length(); s++){
        _dnaMatrix[_indIndex * _nSites + s] = _convert(_str2[s]);
      }
      _indIndex++;
    }
    if(_errCount > 0){
      exit(EXIT_FAILURE);
    }
  }
}

/* Calculate the GH test statistic using counts for current quartet. */
double HyDe::_calcGH(double cp[16][16], int nObs){
  double pxxxx = cp[0][0] + cp[5][5] + cp[10][10] + cp[15][15];
  double pxxxy = cp[0][1] + cp[0][2] + cp[0][3] + cp[5][4]
               + cp[5][6] + cp[5][7] + cp[10][8] + cp[10][9]
               + cp[10][11] + cp[15][12] + cp[15][13] + cp[15][14];
  double pxxyx = cp[0][4] + cp[0][8] + cp[0][12] + cp[5][1]
               + cp[5][9] + cp[5][13] + cp[10][2] + cp[10][6]
               + cp[10][14] + cp[15][3] + cp[15][7] + cp[15][11];
  double pxxyy = cp[0][5] + cp[0][10] + cp[0][15] + cp[5][0]
               + cp[5][10] + cp[5][15] + cp[10][0] + cp[10][5]
               + cp[10][15] + cp[15][0] + cp[15][5] + cp[15][10];
  double pxxyz = cp[0][6] + cp[0][7] + cp[0][9] + cp[0][11]
               + cp[0][13] + cp[0][14] + cp[5][2] + cp[5][3]
               + cp[5][8] + cp[5][11] + cp[5][12] + cp[5][14]
               + cp[10][1] + cp[10][3] + cp[10][4] + cp[10][7]
               + cp[10][12] + cp[10][13] + cp[15][1] + cp[15][2]
               + cp[15][4] + cp[15][6] + cp[15][8] + cp[15][9];
  double pxyxx = cp[1][0] + cp[2][0] + cp[3][0] + cp[4][5]
               + cp[6][5] + cp[7][5] + cp[8][10] + cp[9][10]
               + cp[11][10] + cp[12][15] + cp[13][15] + cp[14][15];
  double pxyxy = cp[1][1] + cp[2][2] + cp[3][3] + cp[4][4]
               + cp[6][6] + cp[7][7] + cp[8][8] + cp[9][9]
               + cp[11][11] + cp[12][12] + cp[13][13] + cp[14][14];
  double pxyxz = cp[1][2] + cp[1][3] + cp[2][1] + cp[2][3]
               + cp[3][1] + cp[3][2] + cp[4][6] + cp[4][7]
               + cp[6][4] + cp[6][7] + cp[7][4] + cp[7][6]
               + cp[8][9] + cp[8][11] + cp[9][8] + cp[9][11]
               + cp[11][8] + cp[11][9] + cp[12][13] + cp[12][14]
               + cp[13][12] + cp[13][14] + cp[14][12] + cp[14][13];
  double pxyyx = cp[1][4] + cp[2][8] + cp[3][12] + cp[4][1]
               + cp[6][9] + cp[7][13] + cp[8][2] + cp[9][6]
               + cp[11][14] + cp[12][3] + cp[13][7] + cp[14][11];
  double pyxxx = cp[4][0] + cp[8][0] + cp[12][0] + cp[1][5]
               + cp[9][5] + cp[13][5] + cp[2][10] + cp[6][10]
               + cp[14][10] + cp[3][15] + cp[7][15] + cp[11][15];
  double pxyyz = cp[1][6] + cp[1][7] + cp[2][9] + cp[2][11]
               + cp[3][13] + cp[3][14] + cp[4][2] + cp[4][3]
               + cp[6][8] + cp[6][11] + cp[7][12] + cp[7][14]
               + cp[8][1] + cp[8][3] + cp[9][4] + cp[9][7]
               + cp[11][12] + cp[11][13] + cp[12][1] + cp[12][2]
               + cp[13][4] + cp[13][6] + cp[14][8] + cp[14][9];
  double pzxyz = cp[8][6] + cp[12][7] + cp[4][9] + cp[12][11]
               + cp[4][13] + cp[8][14] + cp[9][2] + cp[13][3]
               + cp[1][8] + cp[13][11] + cp[1][12] + cp[9][14]
               + cp[6][1] + cp[14][3] + cp[2][4] + cp[14][7]
               + cp[2][12] + cp[6][13] + cp[7][1] + cp[11][2]
               + cp[3][4] + cp[11][6] + cp[3][8] + cp[7][9];
  double pyxzx = cp[4][8] + cp[4][12] + cp[8][4] + cp[8][12]
               + cp[12][4] + cp[12][8] + cp[1][9] + cp[1][13]
               + cp[9][1] + cp[9][13] + cp[13][1] + cp[13][9]
               + cp[2][6] + cp[2][14] + cp[6][2] + cp[6][14]
               + cp[14][2] + cp[14][6] + cp[3][7] + cp[3][11]
               + cp[7][3] + cp[7][11] + cp[11][3] + cp[11][7];
  double pyzxx = cp[6][0] + cp[7][0] + cp[9][0] + cp[11][0]
               + cp[13][0] + cp[14][0] + cp[2][5] + cp[3][5]
               + cp[8][5] + cp[11][5] + cp[12][5] + cp[14][5]
               + cp[1][10] + cp[3][10] + cp[4][10] + cp[7][10]
               + cp[12][10] + cp[13][10] + cp[1][15] + cp[2][15]
               + cp[4][15] + cp[6][15] + cp[8][15] + cp[9][15];
  double pxyzw = cp[1][11] + cp[1][14] + cp[2][7] + cp[2][13]
               + cp[3][6] + cp[3][9] + cp[4][11] + cp[4][14]
               + cp[6][3] + cp[6][12] + cp[7][2] + cp[7][8]
               + cp[8][7] + cp[8][13] + cp[9][3] + cp[9][12]
               + cp[11][1] + cp[11][4] + cp[12][6] + cp[12][9]
               + cp[13][2] + cp[13][8] + cp[14][1] + cp[14][4];

  if(fabs((1.0 / nObs) * (pxxxx + pxxxy + pxxyx + pxxyy + pxxyz + pxyxx + pxyxy + pxyxz
                       +  pxyyx + pyxxx + pxyyz + pzxyz + pyxzx + pyzxx + pxyzw) - 1.0) > 0.005){
    std::cerr << "\n** WARNING: There was a problem counting site patterns ... exiting. **\n" << std::endl;
    //exit(EXIT_FAILURE);
  }

  double p9 = (double) (pxyyx + 0.05) / nObs;
  double p7 = (double) (pxyxy + 0.05) / nObs;
  double p4 = (double) (pxxyy + 0.05) / nObs;
  double obs_invp1 = nObs * (p9 - p7);
  double obs_invp2 = nObs * (p4 - p7);
  if(obs_invp1 == 0){
    obs_invp1 += 1;
    obs_invp2 += 1;
  }
  double obs_var_invp1 = nObs * p9 * (1 - p9) + nObs * p7 * (1 - p7) + 2 * nObs * p9 * p7;
  double obs_var_invp2 = nObs * p4 * (1 - p4) + nObs * p7 * (1 - p7) + 2 * nObs * p4 * p7;
  double obs_cov_invp1_invp2 = -1 * nObs * p9 * p4 + nObs * p9 * p7 + nObs * p7 * p4 + nObs * p7 * (1 - p7);
  double ratio = obs_invp2 / obs_invp1;
  double GH_ts = (obs_invp1) * (ratio) / sqrt(obs_var_invp1 * (pow(ratio, 2.0)) - 2.0
               *  obs_cov_invp1_invp2 * ratio + obs_var_invp2);

  double temp = -99999.9;
  if((p7 > p9) && (p7 > p4)){
    return temp;
  } else if((GH_ts > -99999.9) && (GH_ts < 99999.9)){
    return GH_ts;
  } else {
    return temp;
  }
}

/* Calculate p-value for GH test statistic. */
double HyDe::_calcPvalue(double myZ){
  double a1 = 0.278393;
  double a2 = 0.230389;
  double a3 = 0.000972;
  double a4 = 0.078108;
  double z, erf;
  double myP;

  z = myZ / pow(2,0.5);
  erf = 1.0 - 1.0 / pow(1.0 + a1 * z + a2 * z * z + a3 * z * z * z + a4 * z * z * z * z, 4);
  myP = 2.0 * (1.0 - 0.5 * (1.0 + erf));

  return myP;
}

int HyDe::_getCountMatrix(std::string p1, std::string hyb, std::string p2, double cp[16][16]){
  int nn = 0, resolved = 0;
  for(unsigned i = 0; i < _taxaMap[_outgroup].size(); i++){
    for(unsigned j = 0; j < _taxaMap[p1].size(); j++){
      for(unsigned k = 0; k < _taxaMap[hyb].size(); k++){
        for(unsigned r = 0; r < _taxaMap[p2].size(); r++){
          for(int s = 0; s < _nSites; s++){
            if(_dnaMatrix[_taxaMap[_outgroup][i] * _nSites + s] < 4 &&
               _dnaMatrix[_taxaMap[p1][j] * _nSites + s]  < 4 &&
               _dnaMatrix[_taxaMap[hyb][k] * _nSites + s] < 4 &&
               _dnaMatrix[_taxaMap[p2][r] * _nSites + s]  < 4){
              nn++;
              cp[_dnaMatrix[_taxaMap[_outgroup][i] * _nSites + s] * 4 + _dnaMatrix[_taxaMap[p1][j] * _nSites + s]]
                [_dnaMatrix[_taxaMap[hyb][k] * _nSites + s] * 4 + _dnaMatrix[_taxaMap[p2][r] * _nSites + s]] += 1.0;
            } else {
              resolved = _resolveAmbiguity(_dnaMatrix[_taxaMap[_outgroup][i] * _nSites + s],
                                           _dnaMatrix[_taxaMap[p1][j] * _nSites + s],
                                           _dnaMatrix[_taxaMap[hyb][k] * _nSites + s],
                                           _dnaMatrix[_taxaMap[p2][r] * _nSites + s], cp);
              nn += resolved; /* resolved is 0 or 1. If the site is resolvable it adds to the count, otherwise it doesn't. */
            }
          }

        }
      }
    }
  }
  return nn;
}

int HyDe::_resolveAmbiguity(int out, int p1, int hyb, int p2, double cp[16][16]){
  /* Check to see if any combination of three of the four taxa are ambiguous. */
  double denom = 0.0;
  if((out >= 4 && p1  >= 4 && hyb >= 4) ||
     (out >= 4 && p1  >= 4 && p2  >= 4) ||
     (out >= 4 && hyb >= 4 && p2  >= 4) ||
     (p1  >= 4 && hyb >= 4 && p2  >= 4)){
    return 0;
  } else if(out == 4 || p1 == 4 || hyb == 4 || p2 == 4){ /* Don't allow gaps. */
    return 0;
  } else {
    denom = _baseLookup[out].size() * _baseLookup[p1].size() * _baseLookup[hyb].size() * _baseLookup[p2].size();
    for(unsigned i = 0; i < _baseLookup[out].size(); i++){
      for(unsigned j = 0; j < _baseLookup[p1].size(); j++){
        for(unsigned k = 0; k < _baseLookup[hyb].size(); k++){
          for(unsigned r = 0; r < _baseLookup[p2].size(); r++){
            cp[_baseLookup[out][i] * 4 + _baseLookup[p1][j]]
              [_baseLookup[hyb][k] * 4 + _baseLookup[p2][r]] += 1.0 / denom;
          }
        }
      }
    }
    return 1;
  }
}

void HyDe::_printOut(std::string p1, std::string hyb, std::string p2, double z, double p, double cp[16][16], std::ofstream& out){
  /* These values are now calculated by the program twice. Not sure what the computational cost is. */
  double pxxxx = cp[0][0] + cp[5][5] + cp[10][10] + cp[15][15];
  double pxxxy = cp[0][1] + cp[0][2] + cp[0][3] + cp[5][4]
               + cp[5][6] + cp[5][7] + cp[10][8] + cp[10][9]
               + cp[10][11] + cp[15][12] + cp[15][13] + cp[15][14];
  double pxxyx = cp[0][4] + cp[0][8] + cp[0][12] + cp[5][1]
               + cp[5][9] + cp[5][13] + cp[10][2] + cp[10][6]
               + cp[10][14] + cp[15][3] + cp[15][7] + cp[15][11];
  double pxxyy = cp[0][5] + cp[0][10] + cp[0][15] + cp[5][0]
               + cp[5][10] + cp[5][15] + cp[10][0] + cp[10][5]
               + cp[10][15] + cp[15][0] + cp[15][5] + cp[15][10];
  double pxxyz = cp[0][6] + cp[0][7] + cp[0][9] + cp[0][11]
               + cp[0][13] + cp[0][14] + cp[5][2] + cp[5][3]
               + cp[5][8] + cp[5][11] + cp[5][12] + cp[5][14]
               + cp[10][1] + cp[10][3] + cp[10][4] + cp[10][7]
               + cp[10][12] + cp[10][13] + cp[15][1] + cp[15][2]
               + cp[15][4] + cp[15][6] + cp[15][8] + cp[15][9];
  double pxyxx = cp[1][0] + cp[2][0] + cp[3][0] + cp[4][5]
               + cp[6][5] + cp[7][5] + cp[8][10] + cp[9][10]
               + cp[11][10] + cp[12][15] + cp[13][15] + cp[14][15];
  double pxyxy = cp[1][1] + cp[2][2] + cp[3][3] + cp[4][4]
               + cp[6][6] + cp[7][7] + cp[8][8] + cp[9][9]
               + cp[11][11] + cp[12][12] + cp[13][13] + cp[14][14];
  double pxyxz = cp[1][2] + cp[1][3] + cp[2][1] + cp[2][3]
               + cp[3][1] + cp[3][2] + cp[4][6] + cp[4][7]
               + cp[6][4] + cp[6][7] + cp[7][4] + cp[7][6]
               + cp[8][9] + cp[8][11] + cp[9][8] + cp[9][11]
               + cp[11][8] + cp[11][9] + cp[12][13] + cp[12][14]
               + cp[13][12] + cp[13][14] + cp[14][12] + cp[14][13];
  double pxyyx = cp[1][4] + cp[2][8] + cp[3][12] + cp[4][1]
               + cp[6][9] + cp[7][13] + cp[8][2] + cp[9][6]
               + cp[11][14] + cp[12][3] + cp[13][7] + cp[14][11];
  double pyxxx = cp[4][0] + cp[8][0] + cp[12][0] + cp[1][5]
               + cp[9][5] + cp[13][5] + cp[2][10] + cp[6][10]
               + cp[14][10] + cp[3][15] + cp[7][15] + cp[11][15];
  double pxyyz = cp[1][6] + cp[1][7] + cp[2][9] + cp[2][11]
               + cp[3][13] + cp[3][14] + cp[4][2] + cp[4][3]
               + cp[6][8] + cp[6][11] + cp[7][12] + cp[7][14]
               + cp[8][1] + cp[8][3] + cp[9][4] + cp[9][7]
               + cp[11][12] + cp[11][13] + cp[12][1] + cp[12][2]
               + cp[13][4] + cp[13][6] + cp[14][8] + cp[14][9];
  double pzxyz = cp[8][6] + cp[12][7] + cp[4][9] + cp[12][11]
               + cp[4][13] + cp[8][14] + cp[9][2] + cp[13][3]
               + cp[1][8] + cp[13][11] + cp[1][12] + cp[9][14]
               + cp[6][1] + cp[14][3] + cp[2][4] + cp[14][7]
               + cp[2][12] + cp[6][13] + cp[7][1] + cp[11][2]
               + cp[3][4] + cp[11][6] + cp[3][8] + cp[7][9];
  double pyxzx = cp[4][8] + cp[4][12] + cp[8][4] + cp[8][12]
               + cp[12][4] + cp[12][8] + cp[1][9] + cp[1][13]
               + cp[9][1] + cp[9][13] + cp[13][1] + cp[13][9]
               + cp[2][6] + cp[2][14] + cp[6][2] + cp[6][14]
               + cp[14][2] + cp[14][6] + cp[3][7] + cp[3][11]
               + cp[7][3] + cp[7][11] + cp[11][3] + cp[11][7];
  double pyzxx = cp[6][0] + cp[7][0] + cp[9][0] + cp[11][0]
               + cp[13][0] + cp[14][0] + cp[2][5] + cp[3][5]
               + cp[8][5] + cp[11][5] + cp[12][5] + cp[14][5]
               + cp[1][10] + cp[3][10] + cp[4][10] + cp[7][10]
               + cp[12][10] + cp[13][10] + cp[1][15] + cp[2][15]
               + cp[4][15] + cp[6][15] + cp[8][15] + cp[9][15];
  double pxyzw = cp[1][11] + cp[1][14] + cp[2][7] + cp[2][13]
               + cp[3][6] + cp[3][9] + cp[4][11] + cp[4][14]
               + cp[6][3] + cp[6][12] + cp[7][2] + cp[7][8]
               + cp[8][7] + cp[8][13] + cp[9][3] + cp[9][12]
               + cp[11][1] + cp[11][4] + cp[12][6] + cp[12][9]
               + cp[13][2] + cp[13][8] + cp[14][1] + cp[14][4];
  out << p1 << "\t" << hyb << "\t" << p2 << "\t" << z << "\t" << p << "\t"
      << pxxxx << "\t" << pxxxy << "\t" << pxxyx << "\t" << pxxyy << "\t"
      << pxxyz << "\t" << pxyxx << "\t" << pxyxy << "\t" << pxyxz << "\t"
      << pxyyx << "\t" << pyxxx << "\t" << pxyyz << "\t" << pzxyz << "\t"
      << pyxzx << "\t" << pyzxx << "\t" << pxyzw << std::endl;
}
