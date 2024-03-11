/*
 * This program is licensed granted by STATE UNIVERSITY OF CAMPINAS - UNICAMP ("University")
 * for use of MassCCS software ("the Software") through this website
 * https://github.com/cces-cepid/MassCCS (the "Website").
 *
 * By downloading the Software through the Website, you (the "License") are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about MassCCS please contact: 
 * skaf@unicamp.br (Munir S. Skaf)
 * guido@unicamp.br (Guido Araujo)
 * samuelcm@unicamp.br (Samuel Cajahuaringa)
 * danielzc@unicamp.br (Daniel L. Z. Caetano)
 * zanottol@unicamp.br (Leandro N. Zanotto)
 */

#ifndef MASSCCS_V1_INPUT_H
#define MASSCCS_V1_INPUT_H

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "../rapidjson/document.h"
#include "Constants.h"
#include "omp.h"

using namespace std;

class Input {
private:
  string equipotential_str, gas_buffer_str, short_range_str, long_range_str, long_range, polarizability_str;
  rapidjson::Value atomicParameters;
  rapidjson::Document d;
  
  void parseFile(char const *input);

  void readInputFile(char const *input);

  void checkInput(string option, int type);

public:
  explicit Input(char const *argv);

  void printReadInput();

  unsigned int nProbe;               // numbers of gas buffers
  unsigned int nIter;                // numbers of ccs calculations     
  unsigned int seed;                 // seed       
  unsigned int nthreads;             // numbers of threads, default use all
  string targetFilename;             // molecule target
  string user_ff;                    // user force field
  unsigned int user_ff_flag;         // yes = 1 and not = 0, default is not 
  double dt;                         // time step in fs
  double temperatureTarget;          // temperature in Kelvin
  double skin;                       // skin of linked-cell size
  unsigned int gas_buffer_flag;      // He = 1, N2 = 2, CO2 = 3 and Ar = 4
  unsigned int polarizability_flag;  // yes = 1 and not = 0
  unsigned int equipotential_flag;   // yes = 1 and not = 0
  unsigned int short_range_cutoff;   // yes = 1 and not = 0 for cut lennard-jones interacion
  double lj_cutoff;                  // lennard-jones cutoff   
  unsigned int long_range_flag;      // yes = 1 and not = 0 for apply coulomb interaction
  unsigned int long_range_cutoff;    // yes = 1 and not = 0 for cutoff coulomb interaction
  double coul_cutoff;                // coulomb cutoff      
  double alpha;                      // polarizability 
};

#endif // MASSCCS_V1_INPUT_H
