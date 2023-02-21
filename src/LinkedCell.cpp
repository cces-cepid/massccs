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

#include "headers/LinkedCell.h"

LinkedCell::LinkedCell(MoleculeTarget *moleculeTarget, double a, double b, double c, 
  double lj_cutoff, double skin, unsigned int long_range_flag, unsigned int long_range_cutoff, double coul_cutoff) {
  this->moleculeTarget = moleculeTarget;
  this->a = a;
  this->b = b;
  this->c = c;
  this->lj_cutoff = lj_cutoff;
  this->skin = skin;
  this->long_range_flag = long_range_flag;
  this->long_range_cutoff = long_range_cutoff;
  this->coul_cutoff = coul_cutoff;

  next_neighbor = 0;  
  lx = 2.0*a; 
  ly = 2.0*b; 
  lz = 2.0*c;
  
  calculateNumberOfCells(); 
  calculateAtomsInsideOfCell();
  sortingAtoms();
  calculateCellsNeighbors();  
  print(); 
}

/**
 * Compute the number of cells 
 *
 */
void LinkedCell::calculateNumberOfCells() {
  // number of cell for each axis
  Nx = (int)ceil(lx/(lj_cutoff + skin)); // Nx
  Ny = (int)ceil(ly/(lj_cutoff + skin)); // Ny
  Nz = (int)ceil(lz/(lj_cutoff + skin)); // Nz

  // box simulation domain
  lx = ((double)Nx)*(lj_cutoff + skin); // lx
  ly = ((double)Ny)*(lj_cutoff + skin); // ly
  lz = ((double)Nz)*(lj_cutoff + skin); // lz

  // initialize the LinkedList
  Ncells = Nx*Ny*Nz;
  atoms_inside_cell = new int [Ncells]();
  head_atom_cell = new int[Ncells]();
  atoms_ids = new vector<int> [Ncells]; 
  neighbors1_cells = new int[Ncells]();
  neighbors1_cells_ids = new vector<int> [Ncells];

  if (long_range_flag == 1) {
    next_neighbor = 1;
  }  
  
  if (next_neighbor == 1) {
    neighbors2_cells = new int[Ncells]();
    neighbors2_cells_ids = new vector<int> [Ncells];
  }
  
  for (int i = 0; i < Ncells; i++) {
    atoms_ids[i].push_back(0);
    atoms_ids[i].clear();
    neighbors1_cells_ids[i].push_back(0);
    neighbors1_cells_ids[i].clear();
    if (next_neighbor == 1) {
      neighbors2_cells_ids[i].push_back(0);
      neighbors2_cells_ids[i].clear();
    }  
  } 

  corner.emplace_back(-0.5*lx);
  corner.emplace_back(-0.5*ly);
  corner.emplace_back(-0.5*lz);
}

// calculate the atoms inside each cell
void LinkedCell::calculateAtomsInsideOfCell() {
double pos[3];
int idx, idy, idz; 
int index;
int id;

for (int i = 0; i < moleculeTarget->natoms; i++) {
  pos[0] = moleculeTarget->x[i];
  pos[1] = moleculeTarget->y[i];
  pos[2] = moleculeTarget->z[i];

  calculateIndex(pos,index);

  id = moleculeTarget->id[i];
  atoms_ids[index].push_back(id);
}

for (int i = 0; i < Ncells; i++) { 
  atoms_inside_cell[i] = atoms_ids[i].size(); 
}   

}

// calculate the cell index for specify position
void LinkedCell::calculateIndex(double pos[3], int &index) {
double xi, yi, zi;
int i, j, k;

xi = pos[0] - corner[0];
i = (int)floor(xi/(lj_cutoff + skin));
yi = pos[1] - corner[1];
j = (int)floor(yi/(lj_cutoff + skin));
zi = pos[2] - corner[2];
k = (int)floor(zi/(lj_cutoff + skin));

index = k + Nz*j + Nz*Ny*i;
}

// sorting molecule target atoms
void LinkedCell::sortingAtoms() {
int natoms;
natoms = moleculeTarget->natoms;

string *tmp_atomName;
double *tmp_x;
double *tmp_y;
double *tmp_z;
double *tmp_q;
double *tmp_m;
double *tmp_eps;
double *tmp_sig;

tmp_atomName = new string[natoms]();
tmp_x = new double[natoms]();
tmp_y = new double[natoms]();
tmp_z = new double[natoms]();
tmp_q = new double[natoms]();
tmp_m = new double[natoms]();
tmp_eps = new double[natoms]();
tmp_sig = new double[natoms]();

int iatom;
iatom = 0;
int id;
for (int i = 0; i < Ncells; i++) {
  if (atoms_inside_cell[i]) {
    for (int j = 0; j < atoms_inside_cell[i]; j++) {
       id = atoms_ids[i][j];
       tmp_atomName[iatom] = moleculeTarget->atomName[id];
       tmp_x[iatom] = moleculeTarget->x[id];
       tmp_y[iatom] = moleculeTarget->y[id];
       tmp_z[iatom] = moleculeTarget->z[id];
       tmp_q[iatom] = moleculeTarget->q[id];
       tmp_m[iatom] = moleculeTarget->m[id];
       tmp_eps[iatom] = moleculeTarget->eps[id];
       tmp_sig[iatom] = moleculeTarget->sig[id];
       iatom += 1;  
    }	    
  }	  
}

for (int i = 0; i < natoms; i++) {
   moleculeTarget->atomName[i] = tmp_atomName[i];
   moleculeTarget->x[i] = tmp_x[i];
   moleculeTarget->y[i] = tmp_y[i];
   moleculeTarget->z[i] = tmp_z[i];
   moleculeTarget->q[i] = tmp_q[i];
   moleculeTarget->m[i] = tmp_m[i];
   moleculeTarget->eps[i] = tmp_eps[i];
   moleculeTarget->sig[i] = tmp_sig[i];
}

iatom = 0;
for (int i = 0; i < Ncells; i++) {
   if (atoms_inside_cell[i]) {
     for (int j = 0; j < atoms_inside_cell[i]; j++) {
        atoms_ids[i][j] = iatom;
        iatom += 1;		
     }	     
   } 
}

for (int i = 0; i < Ncells; i++) {
   if (atoms_inside_cell[i]) {
     head_atom_cell[i] = atoms_ids[i][0]; 	   
   } else {
     head_atom_cell[i] = -1; 	   
   }
}

delete [] tmp_atomName;
delete [] tmp_x;
delete [] tmp_y;
delete [] tmp_z;
delete [] tmp_q;
delete [] tmp_m;
delete [] tmp_eps;
delete [] tmp_sig;
}

// calculate the neighbors cell index around specify cell
void LinkedCell::calculateCellsNeighbors() {
// 1: neighbors index calculation
vector<int> indexes(3);
double R1, R2;
double d;
d = lj_cutoff + skin; // cell size
R1 = lj_cutoff/d; // lennard-jones cutoff < 1
if (next_neighbor == 1) R2 = coul_cutoff/d; // coulomb cutoff

vector<vector<int>> indexes_neighbors1{}; 
vector<vector<int>> indexes_neighbors2{};

int ncells, n1, n2;
n1 = (int)ceil(R1); // always 1
ncells = n1;
if (next_neighbor == 1) {
  n2 = (int)ceil(R2); // minimum value is 2 
  ncells = n2;
} 

if (next_neighbor == 1) {
  for (int i = -ncells; i < ncells+1; i++) {  
    for (int j = -ncells; j < ncells+1; j++) {  
      for (int k = -ncells; k < ncells+1; k++) { 
        indexes[0]=i;
        indexes[1]=j;
        indexes[2]=k;     
        if (i==0 && j==0 && k==0) {
          indexes_neighbors1.emplace_back(indexes);           
        } else if (i!=0 && j==0 && k==0) { // axis x
          if (abs(i) <= n1) {
            indexes_neighbors1.emplace_back(indexes);                	     
          } else {
	          indexes_neighbors2.emplace_back(indexes); 	     
          }		     
	      } else if (i==0 && j!=0 && k==0) { // axis y
          if (abs(j) <= n1) {
            indexes_neighbors1.emplace_back(indexes);           
          } else {
            indexes_neighbors2.emplace_back(indexes);
          }          
	      } else if (i==0 && j==0 && k!=0) {	// axis z 
           if (abs(k) <= n1) {
             indexes_neighbors1.emplace_back(indexes);
           } else {
             indexes_neighbors2.emplace_back(indexes);
           }           
        } else if (i==0 && j!=0 && k!=0) { // plane yz
          if ((double)(j*j + k*k) < R2*R2) {
            if ((double)(j*j + k*k) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else if ((double)(pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else {
              indexes_neighbors2.emplace_back(indexes);
	          }
          } else if ((double)(pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R2*R2) {
            if ((double)(j*j + k*k) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else if ((double)(pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else {
              indexes_neighbors2.emplace_back(indexes);
            }  
	        } 
        } else if (i!=0 && j==0 && k!=0) { // plane xz
          if ((double)(i*i + k*k) < R2*R2) {
            if ((double)(i*i + k*k) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else if ((double)(pow(abs(i)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else {
              indexes_neighbors2.emplace_back(indexes);
	          }
          } else if ((double)(pow(abs(i)-1,2) + pow(abs(k)-1,2)) < R2*R2) {
            if ((double)(i*i + k*k) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else if ((double)(pow(abs(i)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else {
              indexes_neighbors2.emplace_back(indexes);
            } 
          }
        } else if (i!=0 && j!=0 && k==0) { // plane xy
          if ((double)(i*i + j*j) < R2*R2) {
            if ((double)(i*i + j*j) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else {
              indexes_neighbors2.emplace_back(indexes);
	          }
          } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2)) < R2*R2) {
            if ((double)(i*i + j*j) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else {
              indexes_neighbors2.emplace_back(indexes);
            }  
          }
        } else {
          if ((double)(i*i + j*j + k*k) < R2*R2) {
            if ((double)(i*i + j*j + k*k) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
	          } else {
              indexes_neighbors2.emplace_back(indexes);
	          }
          } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R2*R2) {
            if ((double)(i*i + j*j + k*k) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
              indexes_neighbors1.emplace_back(indexes);
            } else {
              indexes_neighbors2.emplace_back(indexes);
            }  
          } 
        }   
      }
    }         
  }   
} else {
  for (int i = -ncells; i < ncells+1; i++) {  
    for (int j = -ncells; j < ncells+1; j++) {  
      for (int k = -ncells; k < ncells+1; k++) { 
        indexes[0]=i;
        indexes[1]=j;
        indexes[2]=k;     
        if (i==0 && j==0 && k==0) {
          indexes_neighbors1.emplace_back(indexes);           
        } else if (i!=0 && j==0 && k==0) { // axis x
          if (abs(i) <= n1) {
            indexes_neighbors1.emplace_back(indexes);                	     
          }
	      } else if (i==0 && j!=0 && k==0) { // axis y
          if (abs(j) <= n1) {
            indexes_neighbors1.emplace_back(indexes);           
          }         
	      } else if (i==0 && j==0 && k!=0) {	// axis z 
           if (abs(k) <= n1) {
             indexes_neighbors1.emplace_back(indexes);
           }           
        } else if (i==0 && j!=0 && k!=0) { // plane yz
          if ((double)(j*j + k*k) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
	        } else if ((double)(pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
          } 
        } else if (i!=0 && j==0 && k!=0) { // plane xz
          if ((double)(i*i + k*k) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
	        } else if ((double)(pow(abs(i)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
	        }
        } else if (i!=0 && j!=0 && k==0) { // plane xy
          if ((double)(i*i + j*j) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
	        } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2)) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
	        }
        } else {
          if ((double)(i*i + j*j + k*k) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
	        } else if ((double)(pow(abs(i)-1,2) + pow(abs(j)-1,2) + pow(abs(k)-1,2)) < R1*R1) {
            indexes_neighbors1.emplace_back(indexes);
	        }
        }   
      }
    }         
  }
}

int idx_neighbor;
int idx_cell;

for (int i = 0; i < Nx; i++) {  
  for (int j = 0; j < Ny; j++) {   
    for (int k = 0; k < Nz; k++) {
      idx_cell = k + Nz*j + Nz*Ny*i;     
      // neighbors 1:
      for (int m = 0; m < indexes_neighbors1.size(); m++) {
        indexes[0]= i + indexes_neighbors1.at(m).at(0);
        indexes[1]= j + indexes_neighbors1.at(m).at(1);
        indexes[2]= k + indexes_neighbors1.at(m).at(2);
        if (-1 < indexes[0] && indexes[0] < Nx && -1 < indexes[1] && indexes[1] < Ny && -1 < indexes[2] && indexes[2] < Nz) {
	        idx_neighbor = indexes[2] + Nz*indexes[1] + Nz*Ny*indexes[0]; 	    
          if (atoms_inside_cell[idx_neighbor] != 0) neighbors1_cells_ids[idx_cell].emplace_back(idx_neighbor);
        }		    
      }
      if (next_neighbor == 1) {
        // neighbors 2:
        for (int m = 0; m < indexes_neighbors2.size(); m++) {
          indexes[0]= i + indexes_neighbors2.at(m).at(0);
          indexes[1]= j + indexes_neighbors2.at(m).at(1);
          indexes[2]= k + indexes_neighbors2.at(m).at(2);
          if (-1 < indexes[0] && indexes[0] < Nx && -1 < indexes[1] && indexes[1] < Ny && -1 < indexes[2] && indexes[2] < Nz) {
            idx_neighbor = indexes[2] + Nz*indexes[1] + Nz*Ny*indexes[0];	    
            if (atoms_inside_cell[idx_neighbor] != 0) neighbors2_cells_ids[idx_cell].emplace_back(idx_neighbor);
          }
        }
      }   
    }
  }
}   

for (int i = 0; i < Ncells; i++) {  // <- loop: i < Nx
  neighbors1_cells[i] = neighbors1_cells_ids[i].size();     
  if (next_neighbor == 1) neighbors2_cells[i] = neighbors2_cells_ids[i].size();
}

}

// print information
void LinkedCell::print() {
int filled_cells, empty_cells, average_atoms_cells, maximum_atoms_cells, minimum_atoms_cells;

average_atoms_cells = 0;
maximum_atoms_cells = 0;
minimum_atoms_cells = moleculeTarget->natoms;
empty_cells = 0;
filled_cells = 0;
for (int i = 0; i < Ncells; i++) {  
   if (atoms_inside_cell[i] == 0) {
     empty_cells += 1;
   } else {
     average_atoms_cells += atoms_inside_cell[i];	
     filled_cells += 1;
     if (atoms_inside_cell[i] > maximum_atoms_cells) {
       maximum_atoms_cells = atoms_inside_cell[i]; 		   
     }		   
     if (atoms_inside_cell[i] < minimum_atoms_cells) {
       minimum_atoms_cells = atoms_inside_cell[i];
     }   	   
   }	 
}

cout << "*********************************************************"
            << endl;
cout << "Linked-cell: " << endl;
cout << "Numbers of cells: " << Ncells << endl;
cout << "Nx: " << Nx << " Ny: " << Ny << " Nz: " << Nz << endl;
cout << "Filled cells: " << filled_cells << endl;
cout << "Empty cells: " << empty_cells << endl;
cout << "Average atoms per cell: " << ((float)average_atoms_cells)/((float)filled_cells) << endl;
cout << "Maximum atoms per cell: " << maximum_atoms_cells << endl;
cout << "Minimum atoms per cell: " << minimum_atoms_cells << endl;
cout << "Simulation box: " << endl;
cout << "lx: " << lx << " Ang" << endl;
cout << "ly: " << ly << " Ang" << endl;
cout << "lz: " << lz << " Ang" << endl;
}
