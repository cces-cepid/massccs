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

#ifndef MASSCCS_V1_MOLECULETARGET_H
#define MASSCCS_V1_MOLECULETARGET_H

#include "MoleculeTarget.h"
#include "GasBuffer.h"
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>
#include "Math.h"
#include <algorithm>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <iostream>

using namespace std;

class MoleculeTarget {
private:
  string filename;
  string extension;
  string user_ff;
  unsigned int gas_buffer_flag;
  unsigned int user_ff_flag;
  unsigned int force_type;
  void readXYZfile(string &filename);
  void readPQRfile(string &filename);
  void readUserFF(string &user_ff);
  void defaultFF();
  void calculateCenterOfMass(double(&)[3]);
  void moveToCenterOfMass(double[3]);
  void calculateMoleculeRadius();
  void orientateMolecule();
  void setInertia();
  void getEigenvalues(double[3][3], double(&)[3]);
  void getEigenvectors(double[3][3], double[3], double(&)[3][3]);
  double delta(int, int);
  double inertia[3][3];
  double inertiaValues[3];
  double inertiaVectors[3][3];
  vector<double> assignedParameter(string chemical);
  void printFF();
  unsigned int nparameters;
  string *user_atomName;
  double *user_m;
  double *user_eps;
  double *user_sig;
  bool diagonal; 

  void print();

public:
  MoleculeTarget(string &filename, unsigned int gas_buffer_flag, string &user_ff, unsigned int user_ff_flag, unsigned int force_type);

  unsigned int natoms;

  int *id;
  string *atomName;
  double *x;
  double *y;
  double *z;
  double *q;
  double *m;
  double *eps;
  double *sig;

  double *eps_central;
  double *sig_central;

  double moleculeRadius = 0.0, mass = 0.0, Q = 0.0;
  double rcm[3];

  double maxX, maxY, maxZ;
};

#endif // MASSCCS_V1_MOLECULETARGET_H
