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

#ifndef MASSCCS_V1_EQUIPOTENTIAL_H
#define MASSCCS_V1_EQUIPOTENTIAL_H

#include "MoleculeTarget.h"
#include "GasBuffer.h"
#include "omp.h"
#include <vector>
#include "Math.h"
#include <algorithm>
#include "Constants.h"

using namespace std;

class Equipotential {

private:
  vector<vector<double>> boundryPoints{};
  MoleculeTarget *moleculeTarget;
  GasBuffer *gas;
  unsigned int polarizability_flag;
  unsigned int gas_buffer_flag;
  double epsilon_gas, sigma_gas;
  double aUncert, bUncert, cUncert, rmsd;
  double cov_aa, cov_bb, cov_cc, cov_ab, cov_ac, cov_bc;
  double avgAxes, avgAxesUncert;
  double surfArea, surfAreaUncert;
  double minPhi, minCosT, minSinT;
  double aMin, bMin, cMin;
  double aUncertMin, bUncertMin, cUncertMin;
  double enlargeAmount;
  double v_mb, v2_mb, sig_mb, v_min_mb, Ek_min;
  double temperature;
  double mu;
  double alpha;
  vector<vector<double>> minmax{};

public:
  Equipotential(MoleculeTarget *moleculeTarget, GasBuffer *gas, unsigned int long_range_flag, double temperature, double mu, double alpha, unsigned int gas_buffer_flag);
  ~Equipotential();

  double maxX, maxY, maxZ;
  double a, b, c; // axes of ellipsoid

  void equipotentialpoints();

  void ellipsoid();

  void print();  

  //double potential(vector<double> pos); 
  double potential_He(vector<double> pos);
  double potential_N2(vector<double> pos);
  double potential_CO2(vector<double> pos);

  void enlargeEllipsoidBoundry();
};

#endif // MASSCCS_V1_EQUIPOTENTIAL_H
