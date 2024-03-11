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

#ifndef MASSCCS_V1_FORCE_H
#define MASSCCS_V1_FORCE_H

#include "MoleculeTarget.h"
#include "GasBuffer.h"
#include <cmath>
#include <vector>
#include "LinkedCell.h"
#include "Constants.h"

class Force {
private:
  double lj_cutoff;
  double alpha;
  double coul_cutoff;
  LinkedCell *linkedcell;
  int Nx, Ny, Nz;
  MoleculeTarget *moleculeTarget;
  
public:
 
  Force(MoleculeTarget *moleculeTarget, LinkedCell *linkedcell, double lj_cutoff, double alpha, double coul_cutoff); 

  ~Force();

  void lennardjones(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_coulomb(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_coulomb_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void coulomb(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void coulomb_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_induced_dipole(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_induced_dipole_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void coulomb_induced_dipole(GasBuffer *gas, int iatom, vector<double> n, vector<double> &f, double &Up);
  void coulomb_induced_dipole_LC(GasBuffer *gas, int iatom, vector<double> n, vector<double> &f, double &Up);
  void coulomb_induced_dipole_iso(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void coulomb_induced_dipole_iso_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void nitrogen2_LC(GasBuffer *gas, vector<double> &f, double &Up);
  void nitrogen(GasBuffer *gas, vector<double> &f, double &Up);
  // for carbon of CO2 molecule
  void lennardjones_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_LC_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_coulomb_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_coulomb_LC_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_coulomb_induced_dipole_iso_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up);
  void lennardjones_coulomb_induced_dipole_iso_LC_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up);

  //double Switch(double rIn, double rOut, double r);
  //double DSwitch(double rIn, double rOut, double r);
};

#endif 
