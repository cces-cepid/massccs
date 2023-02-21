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

#ifndef MASSCCS_V1_LINKEDCELL_H
#define MASSCCS_V1_LINKEDCELL_H

#include "MoleculeTarget.h"
#include <cmath>
#include <vector>

using namespace std;

class LinkedCell {
private:
  vector<double> corner{};
  double lj_cutoff;
  unsigned int long_range_flag;
  unsigned int long_range_cutoff;
  int next_neighbor;
  double coul_cutoff;
  double skin;
  MoleculeTarget *moleculeTarget;

  void calculateNumberOfCells();
  void calculateAtomsInsideOfCell();
  void sortingAtoms();
  void calculateCellsNeighbors();
  void print();

public:
  LinkedCell(MoleculeTarget *moleculeTarget, double a, double b, double c,
	  double lj_cutoff, double skin, unsigned int long_range_flag, unsigned int long_range_cutoff, double coul_cutoff);

  int Nx, Ny, Nz, Ncells;
  double a, b, c;
  int *atoms_inside_cell;
  int *head_atom_cell;
  vector<int> *atoms_ids;
  int *neighbors1_cells; 
  vector<int> *neighbors1_cells_ids;
  int *neighbors2_cells;
  vector<int> *neighbors2_cells_ids;
  double lx, ly, lz;
  void calculateIndex(double [3], int &index);
};

#endif // MASSCCS_V1_LINKEDCELL_H
