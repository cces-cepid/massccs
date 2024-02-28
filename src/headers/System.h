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

#ifndef MASSCCS_V1_SYSTEM_H
#define MASSCCS_V1_SYSTEM_H

#include <fstream>
#include <map>

#include "Input.h"
#include "Time.h"
#include "omp.h"
#include <algorithm>
#include <vector>
#include "RandomNumber.h"
#include "Math.h"
#include "MoleculeTarget.h"
#include "GasBuffer.h"
#include "Equipotential.h"
#include "LinkedCell.h"
#include "Force.h"
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <iostream>

using namespace std;

class System {
private:
  unsigned int seed, nProbe, nIter, equipotential_flag, gas_buffer_flag, nthreads;
  unsigned int short_range_cutoff, long_range_flag, long_range_cutoff, polarizability_flag, user_ff_flag;
  double temperatureTarget, dt, skin;
  double lj_cutoff;
  double coul_cutoff;
  unsigned int force_type; 
  string targetFilename, user_ff;
  Input *input;
  RandomNumber *mt{};
  MoleculeTarget *moleculeTarget;
  GasBuffer *gas;
  Equipotential *equipotential;
  LinkedCell *linkedcell;

  double ccs{}, ccs2{};
  double CCS_ave, CCS_err;
  double mu; 
  double alpha;
  double Inertia;
  double d_bond;

  double bmax;
  double a, b, c;
  double lx, ly, lz;

  void setup(GasBuffer *gasProbe,bool &hit, double rndVal1, double rndVal2, double rndVal3, double rndVal4, double rndVal5, 
double rndVal6, double rndVal7, double rndVal8, double rndVal9);
  double velDistr(double v, double m, double temperature);
  double velGenerator(double m, double temperature, double sd);
  double KineticEnergy(double m, vector<double> v);
  double anglevec(vector<double> vi, vector<double> vf);
  void rotate(vector<double> &r, vector<double> &v, vector<double> angles);
  void rotate_gas(vector<double> &v, double theta, double phi);
  void geometric_ellipsoid();

  void first_half_verlet_constrained(vector<double> &ri, vector<double> &rj, vector<double> &vi, vector<double> &vj,
 vector<double> fi, vector<double> fj, double mi, double mj, double dt, double d2ij);

  void second_half_verlet_constrained(vector<double> ri, vector<double> rj, vector<double> &vi, vector<double> &vj,
 vector<double> fi, vector<double> fj, double mi, double mj, double dt);  

public:

  explicit System(char *inputFilename);
  
  void run_He(GasBuffer *gasProbe, bool &success, double &chi, double dt, Force *force);
  void run_N2(GasBuffer *gasProbe, bool &success, double &chi, double dt, Force *force);
  void run_CO2(GasBuffer *gasProbe, bool &success, double &chi, double dt, Force *force);
  //void run_N2_one_site(GasBuffer *gasProbe, bool &success, double &chi, double dt, Force *force);

  ~System();
};

#endif 
