/*
 * This program is licensed granted by STATE UNIVERSITY OF CAMPINAS - UNICAMP ("University")
 * for use of HIGH PERFORMANCE COLLISION CROSS SECTION LINKED-CELL - HPCCS-LC ("the Software") through this website
 * --------------------------------------- (the "Website").
 *
 * By downloading the Software through the Website, you (the "Licensee") are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about HPCCS-LC please contact: 
 * skaf@iqm.unicamp.br (Munir Skaf)
 * samuelcajahuaringa@gmail.com (Samuel Cajahuaringa)
 * leandro.zanotto@gmail.com (Leandro Zanotto)
 * daniel.lzcaetano@gmail.com (Daniel Caetano)
 */

#include "headers/Force.h"

Force::Force(MoleculeTarget *moleculeTarget, LinkedCell *linkedcell, double lj_cutoff, double alpha, double coul_cutoff) {
this->moleculeTarget = moleculeTarget; 
this->linkedcell = linkedcell;           
this->lj_cutoff = lj_cutoff;                     
this->alpha = alpha;
this->coul_cutoff = coul_cutoff;
}

Force::~Force(){	
}

/*
 * Compute the lennard jones force and potential using linked-cell
 */
void Force::lennardjones_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy, fz, U, Ulj, flj, Ulj_cut;    
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, r;
double r2inv, r6inv, rc6inv;
double lj1, lj2, lj3, lj4;
double s1, s2;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;
    
int index; 
linkedcell->calculateIndex(r_probe,index); 
  
if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);	    
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
} 	    

int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;  
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];  
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);

    if (r < lj_cutoff) {
      epsilon = moleculeTarget->eps[target_id]; 
      sigma = moleculeTarget->sig[target_id]; 
  
      r2inv = 1.0/r2;
      r6inv = r2inv*r2inv*r2inv;
      lj1 = 4.0*epsilon*pow(sigma,6.0);
      lj2 = lj1*pow(sigma,6.0);
      Ulj = r6inv*(lj2*r6inv - lj1);
      rc6inv = 1.0/pow(lj_cutoff,6);    
      Ulj_cut = rc6inv*(lj2*rc6inv - lj1); 
      lj3 = 6.0*lj1;
      lj4 = 12.0*lj2;
      flj = r6inv*(lj4*r6inv - lj3)*r2inv;

      U += Ulj - Ulj_cut;
      fx += flj*dx;
      fy += flj*dy;
      fz += flj*dz;
    }
  }	    
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;
return;
}

/*
 * Compute the lennard jones force and potential
 */

void Force::lennardjones(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2,r;
double r2inv,r6inv;
double lj1,lj2,lj3,lj4;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];

Ulj = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);

  //if (r < lj_cutoff) {
  epsilon = moleculeTarget->eps[i];
  sigma = moleculeTarget->sig[i];
     
  r2inv = 1.0/r2;
  r6inv = r2inv*r2inv*r2inv;
  lj1 = 4.0*epsilon*pow(sigma,6.0);
  lj2 = lj1*pow(sigma,6.0);
  Ulj += r6inv*(lj2*r6inv - lj1);

  lj3 = 6.0*lj1;
  lj4 = 12.0*lj2;
  flj = r6inv*(lj4*r6inv - lj3)*r2inv;
  fx += flj*dx;
  fy += flj*dy;
  fz += flj*dz;
  //}
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = Ulj;

return;
}

/**
 * Compute the lennard jones and coulomb interactions
 **/
void Force::lennardjones_coulomb(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, r;
double rinv, r2inv, r6inv;
double lj1, lj2, lj3, lj4;
double qi, qj;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

double Ucoul, fcoul, s1, s2;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  //if (r < lj_cutoff) {
  epsilon = moleculeTarget->eps[i];
  sigma = moleculeTarget->sig[i];
     
  r2inv = 1.0/r2; 
  r6inv = r2inv*r2inv*r2inv;
  lj1 = 4.0*epsilon*pow(sigma,6.0);
  lj2 = lj1*pow(sigma,6.0);
  Ulj = r6inv*(lj2*r6inv - lj1);
  lj3 = 6.0*lj1;
  lj4 = 12.0*lj2;
  flj = r6inv*(lj4*r6inv - lj3)*r2inv;
  
  U += Ulj;
  fx += flj*dx;
  fy += flj*dy;
  fz += flj*dz;
  //}

  //if (r < coul_cutoff) {
  qj = moleculeTarget->q[i];
  rinv = 1.0/r;
  Ucoul = qi*qj*rinv*KCOUL;
  r2inv = 1.0/r2;
  fcoul = Ucoul*r2inv; 
  U += Ucoul;
  fx += fcoul*dx;
  fy += fcoul*dy;
  fz += fcoul*dz; 
  //}
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;

return;
}

/**
 * Compute the coulomb interactions
 **/
void Force::coulomb(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz;
double dx, dy, dz;
double r2,r;
double rinv,r2inv;
double Ucoul, fcoul, U;
double qi, qj;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

fx = 0.0;
fy = 0.0;
fz = 0.0;
U = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);
  
  //if(r < coul_cutoff) {
  qj = moleculeTarget->q[i];
  rinv = 1.0/r;
  Ucoul = qi*qj*rinv*KCOUL;
  r2inv = 1.0/r2;
  fcoul = Ucoul*r2inv;

  U += Ucoul;
  fx += fcoul*dx;
  fy += fcoul*dy;
  fz += fcoul*dz;
  //}
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;

return;
}

/**
 * Compute the lennard jones and induced dipole interactions using linked-cell (apply for Helium)
 **/
void Force::lennardjones_induced_dipole_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U, Ulj_cut;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, x2, y2, z2, r;
double r2inv, r6inv, rc6inv;
double lj1, lj2, lj3, lj4;
double q, qr3inv, qr5inv, qrc;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double r3inv, r5inv, r7inv, r9inv;
double rc3inv, smooth_factor;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

int index;
linkedcell->calculateIndex(r_probe,index);

if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
}

Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;

double s1, s2, Exi, Eyi, Ezi, Exxi, Eyyi, Ezzi, Exyi, Exzi, Eyzi;

// calculation lennard-jones and induced dipole interactions on the first neighbors cells
int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);
    r2inv = 1.0/r2;
    
    // lennard-jones interaction
    if (r < lj_cutoff) {
     epsilon = moleculeTarget->eps[target_id];
     sigma = moleculeTarget->sig[target_id];
     r6inv = r2inv*r2inv*r2inv;
     lj1 = 4.0*epsilon*pow(sigma,6.0);
     lj2 = lj1*pow(sigma,6.0);
     Ulj = r6inv*(lj2*r6inv - lj1);
     rc6inv = 1.0/pow(lj_cutoff,6);
     Ulj_cut = rc6inv*(lj2*rc6inv - lj1);

     lj3 = 6.0*lj1;
     lj4 = 12.0*lj2;
     flj = r6inv*(lj4*r6inv - lj3)*r2inv;
       
     U += Ulj - Ulj_cut;
     fx += flj*dx;
     fy += flj*dy;
     fz += flj*dz;
    }

    // ion-induced dipole interaction
    if (r < coul_cutoff) {
     r3inv = 1.0/r*r2inv;
     r5inv = r3inv*r2inv;
     q = moleculeTarget->q[target_id];
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = q*r3inv*smooth_factor;
     qr5inv = -3.0*q*r5inv*smooth_factor;
     qrc = -3.0*q*rc3inv*r2inv; 

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;
    }
  }
}

// calculation induced dipole interaction on the second neighbors cells
neighborscells = linkedcell->neighbors2_cells[index];
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors2_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);
    r2inv = 1.0/r2;

    // ion-induced dipole interaction
    if (r < coul_cutoff) {
     r3inv = 1.0/r*r2inv;
     r5inv = r3inv*r2inv;
     q = moleculeTarget->q[target_id];
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = q*r3inv*smooth_factor;
     qr5inv = -3.0*q*r5inv*smooth_factor;
     qrc = -3.0*q*rc3inv*r2inv;

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;
    }
  }
}

f[0] = fx + alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
f[1] = fy + alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
f[2] = fz + alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

Up = U - 0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
return;
}

/*
 * Compute the lennard jones and induced dipole interactions (Hellium atom)
 */

void Force::lennardjones_induced_dipole(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, x2, y2, z2, r;
double r2inv, r6inv;
double lj1, lj2, lj3, lj4;
double q, qr3inv, qr5inv;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double r3inv, r5inv, r7inv, r9inv;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];

Ulj = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  epsilon = moleculeTarget->eps[i];
  sigma = moleculeTarget->sig[i];

  r2inv = 1.0/r2;
  r6inv = r2inv*r2inv*r2inv;
  lj1 = 4.0*epsilon*pow(sigma,6.0);
  lj2 = lj1*pow(sigma,6.0);
  Ulj += r6inv*(lj2*r6inv - lj1);
  
  lj3 = 6.0*lj1;
  lj4 = 12.0*lj2;
  flj = r6inv*(lj4*r6inv - lj3)*r2inv;
  fx += flj*dx;
  fy += flj*dy;
  fz += flj*dz;

  q =  moleculeTarget->q[i];
  r3inv = 1.0/r*r2inv;
  r5inv = r3inv*r2inv;
  qr3inv = q*r3inv;
  qr5inv = q*r5inv;

  Ex += dx * qr3inv;
  Ey += dy * qr3inv;
  Ez += dz * qr3inv;

  Exx += qr3inv - 3.0*dx*dx*qr5inv;
  Eyy += qr3inv - 3.0*dy*dy*qr5inv;
  Ezz += qr3inv - 3.0*dz*dz*qr5inv;

  Exy += - 3.0 * dx * dy * qr5inv;
  Exz += - 3.0 * dx * dz * qr5inv;
  Eyz += - 3.0 * dy * dz * qr5inv;
}

f[0] = fx + alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
f[1] = fy + alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
f[2] = fz + alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

Up = Ulj - 0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
return;
}

/*double Force::Switch(double rIn, double rOut, double r) {
double a1, a2, a3;

a1 = pow(rOut*rOut -r*r,2);
a2 = rOut*rOut + 2.0*r*r - 3.0*rIn*rIn;
a3 = pow(rOut*rOut - rIn*rIn,3);

return a1*a2/a3;
}

double Force::DSwitch(double rIn, double rOut, double r) {
double a1, a2, a3;

a1 = rOut*rOut -r*r;
a2 = r*r - rIn*rIn;
a3 = pow(rOut*rOut - rIn*rIn,3);
return -12.0*a1*a2/a3;
}*/

/*
 * Compute the lennard jones and coulomb interactions using linked-cell (<- N of N2 molecule)
 */
void Force::lennardjones_coulomb_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2,r;
double rinv,r2inv,r6inv;
double lj1,lj2,lj3,lj4;
double qi, qj;
double Ucoul, fcoul, s1, s2;
int index;
double Ulj_cut, rc6inv, Ucoul_shift;


r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

linkedcell->calculateIndex(r_probe,index);

if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
}

int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);

    if (r < lj_cutoff) {
     epsilon = moleculeTarget->eps[target_id];
     sigma = moleculeTarget->sig[target_id]; 
     r2inv = 1.0/r2;
     r6inv = r2inv*r2inv*r2inv;
     lj1 = 4.0*epsilon*pow(sigma,6.0);
     lj2 = lj1*pow(sigma,6.0);
     Ulj = r6inv*(lj2*r6inv - lj1);
     rc6inv = 1.0/pow(lj_cutoff,6);
     Ulj_cut = rc6inv*(lj2*rc6inv - lj1);
     lj3 = 6.0*lj1;
     lj4 = 12.0*lj2;
     flj = r6inv*(lj4*r6inv - lj3)*r2inv;		

     U += Ulj - Ulj_cut;
     fx += flj*dx;
     fy += flj*dy;
     fz += flj*dz;
    }

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));
     	
     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz; 
    } 	      
  }
}

// calculation induced dipole interaction on the second neighbors cells
neighborscells = linkedcell->neighbors2_cells[index];
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors2_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz;
    }
  }
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;

return;
}

/**
 * Compute coulomb interaction using linked-cell (<- central charge of N2 molecule)
 **/
void Force::coulomb_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy, fz, Ulj, flj, U;
double dx, dy, dz;
double r2, r;
double rinv, r2inv, r6inv;
double lj1, lj2, lj3, lj4;
double qi, qj;
double Ucoul_shift;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

double Ucoul, fcoul, s1, s2;

int index;
linkedcell->calculateIndex(r_probe,index);

if (index >= linkedcell->Ncells && index < 0) {
  printf("outside cell simulation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
}

int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r = sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz; 
    } 	      
  }
}

// calculation induced dipole interaction on the second neighbors cells
neighborscells = linkedcell->neighbors2_cells[index];
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors2_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r = sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz;
    }
  }
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;

return;
}

/*
 * Compute the coulomb and induced dipole interactions with isotropy polarizability (<- central charge of N2 molecule)
 */

void Force::coulomb_induced_dipole_iso(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, U;
double dx, dy, dz;
double r2, x2, y2, z2, r;
double r2inv, r6inv;
double fcoul, Ucoul;
double qi, qj, qr3inv, qr5inv, qrc;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double rinv, r3inv, r5inv, r7inv, r9inv, rc3inv;
double U_ind, f_ind[3];

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  qj = moleculeTarget->q[i];
  rinv = 1.0/r;
  Ucoul = qi*qj*rinv*KCOUL;
  r2inv = 1.0/r2;
  fcoul = Ucoul*r2inv;

  U += Ucoul;
  fx += fcoul*dx;
  fy += fcoul*dy;
  fz += fcoul*dz;

  r3inv = rinv*r2inv;
  r5inv = r3inv*r2inv;
  qr3inv = qj*r3inv;
  qr5inv = qj*r5inv;

  Ex += dx * qr3inv;
  Ey += dy * qr3inv;
  Ez += dz * qr3inv;

  Exx += qr3inv - 3.0*dx*dx*qr5inv;
  Eyy += qr3inv - 3.0*dy*dy*qr5inv;
  Ezz += qr3inv - 3.0*dz*dz*qr5inv;

  Exy += - 3.0 * dx * dy * qr5inv;
  Exz += - 3.0 * dx * dz * qr5inv;
  Eyz += - 3.0 * dy * dz * qr5inv;
}

f[0] = fx + alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
f[1] = fy + alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
f[2] = fz + alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

Up = U - 0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
return;
}

/*
 * Compute the coulomb and induced dipole interactions with anisotropy polarizability (<- central charge of N2 molecule)
 */

void Force::coulomb_induced_dipole_aniso(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, U;
double dx, dy, dz;
double r2, x2, y2, z2, r;
double r2inv, r6inv;
double fcoul, Ucoul;
double qi, qj, qr3inv, qr5inv, qrc;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double rinv, r3inv, r5inv, r7inv, r9inv, rc3inv;
double U_ind, f_ind[3];
double EE[3][3], E[3], alpha_tensor[3][3];
double alpha_parallel, alpha_perpendicular;
double un[3];
double theta, phi;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

// orientation of nitrogen molecule
un[0] = gas->x[0] - gas->x[1];
un[1] = gas->y[0] - gas->y[1];
un[2] = gas->z[0] - gas->z[1];

double un2;
un2 = gas->d;

un[0] /= un2;
un[1] /= un2;
un[2] /= un2;

theta = acos(un[2]);
phi = atan2(un[1],un[0]);

//polarizability (axial)  = 2.19609742e-30 in m3
//polarizability (radial) = 1.51148405e-30 in m3
alpha_parallel = 2.19609742; 
alpha_perpendicular = 1.51148405;
alpha_parallel *= ALPHA_TO_KCAL_MOL;
alpha_perpendicular *= ALPHA_TO_KCAL_MOL;

alpha_tensor[0][0] = alpha_perpendicular*pow(sin(phi),2) + pow(cos(phi),2)*(alpha_perpendicular*pow(cos(theta),2)+alpha_parallel*pow(sin(theta),2));
alpha_tensor[1][1] = alpha_perpendicular*pow(cos(phi),2) + pow(sin(phi),2)*(alpha_perpendicular*pow(cos(theta),2)+alpha_parallel*pow(sin(theta),2));
alpha_tensor[2][2] = alpha_parallel*pow(cos(theta),2) + alpha_perpendicular*pow(sin(theta),2);
alpha_tensor[0][1] = (alpha_parallel-alpha_perpendicular)*cos(phi)*sin(phi)*pow(sin(theta),2);
alpha_tensor[1][0] = alpha_tensor[0][1];
alpha_tensor[0][2] = (alpha_parallel-alpha_perpendicular)*cos(phi)*cos(theta)*sin(theta);
alpha_tensor[2][0] = alpha_tensor[0][2];
alpha_tensor[1][2] = (alpha_parallel-alpha_perpendicular)*cos(theta)*sin(phi)*sin(theta);
alpha_tensor[2][1] = alpha_tensor[1][2];

/*for (int i = 0; i < 3; i++) {
  for (int j = 0; j < 3; j++) {
    alpha_d[i][j] = 0.0;
    alpha_tensor[i][j] = 0.0;
  }
}

alpha_d[0][0] = alpha_radial;
alpha_d[1][1] = alpha_radial;
alpha_d[2][2] = alpha_axial;


R[0][0] = cos(phi)*cos(theta);
R[1][0] = sin(phi)*cos(theta);
R[2][0] = -sin(theta);

R[0][1] = -sin(phi);
R[1][1] = cos(phi);
R[2][1] = 0.0;

R[0][2] = cos(phi)*sin(theta);
R[1][2] = sin(phi)*sin(theta);
R[2][2] = cos(theta);

for (int i = 0; i < 3; i++) {
  for (int j = 0; j < 3; j++) {
    RT[i][j] = R[j][i];
  }
}

for (int i = 0; i < 3; i++) {
  for (int j = 0; j < 3; j++) {
    for (int k = 0; k < 3; k++) {
      for (int m = 0; m < 3; m++) {
        alpha_tensor[i][j] += R[i][k]*alpha_d[k][m]*RT[m][j];
      }
    }
  }
}*/

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  qj = moleculeTarget->q[i];
  rinv = 1.0/r;
  Ucoul = qi*qj*rinv*KCOUL;
  r2inv = 1.0/r2;
  fcoul = Ucoul*r2inv;

  U += Ucoul;
  fx += fcoul*dx;
  fy += fcoul*dy;
  fz += fcoul*dz;

  r3inv = rinv*r2inv;
  r5inv = r3inv*r2inv;
  qr3inv = qj*r3inv;
  qr5inv = qj*r5inv;

  Ex += dx * qr3inv;
  Ey += dy * qr3inv;
  Ez += dz * qr3inv;

  Exx += qr3inv - 3.0*dx*dx*qr5inv;
  Eyy += qr3inv - 3.0*dy*dy*qr5inv;
  Ezz += qr3inv - 3.0*dz*dz*qr5inv;

  Exy += - 3.0 * dx * dy * qr5inv;
  Exz += - 3.0 * dx * dz * qr5inv;
  Eyz += - 3.0 * dy * dz * qr5inv;
}

// induced-dipole contribution
E[0] = Ex;
E[1] = Ey;
E[2] = Ez;

EE[0][0] = Exx; //Exx
EE[1][1] = Eyy; //Eyy
EE[2][2] = Ezz; //Ezz

EE[0][1] = Exy; //Exy
EE[1][0] = EE[0][1]; //Eyx
EE[0][2] = Exz; //Exz
EE[2][0] = EE[0][2]; //Ezx
EE[1][2] = Eyz; //Eyz
EE[2][1] = EE[1][2]; //Ezy

U_ind = 0.0;
for (int i = 0; i < 3; i++) {
  for (int j = 0; j < 3; j++) {
    U_ind += -0.5*alpha_tensor[i][j]*E[j]*E[i];
  }
}

Up = U + U_ind;

for (int k = 0; k < 3; k++) {
  f_ind[k] = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      f_ind[k] += 0.5*alpha_tensor[i][j]*(EE[j][k]*E[i] + E[j]*EE[i][k]);
    }
  }
}

f[0] = fx + f_ind[0];
f[1] = fy + f_ind[1];
f[2] = fz + f_ind[2];
/*
double Uiso, fiso[3];

fiso[0] = alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
fiso[1] = alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
fiso[2] = alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

Uiso = - 0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);

//cout << "ratio aniso/iso: " << U_ind/Uiso << "  " << U_ind << "  "  << Uiso << endl;
theta *= RAD_TO_DEGREE;
phi *= RAD_TO_DEGREE;
if (phi < 0.0) phi += 360.;

double f1 = sqrt(f_ind[0]*f_ind[0]+f_ind[1]*f_ind[1]+f_ind[2]*f_ind[2]);
double f2 = sqrt(fiso[0]*fiso[0]+fiso[1]*fiso[1]+fiso[2]*fiso[2]);
double beta = acos((f_ind[0]*fiso[0]+f_ind[1]*fiso[1]+f_ind[2]*fiso[2])/(f1*f2))*RAD_TO_DEGREE;

cout<< theta << "  " << phi << "  " << f1 << "  " << f2 << "  " << beta << endl;//f_ind[0]/fiso[0] << "  " << f_ind[1]/fiso[1] << "  " << f_ind[2]/fiso[2] << "  " << endl;
*/
return;
}

/*
 * Compute the coulomb and induced dipole interactions with isotropy polarizability using linked-cell list (<- central charge of N2 molecule)
 */

void Force::coulomb_induced_dipole_iso_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, U;
double dx, dy, dz;
double r2, x2, y2, z2, r;
double r2inv, r6inv, rinv;
double qi, qj, qr3inv, qr5inv, qrc;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double Exi, Eyi, Ezi, Exxi, Exyi, Exzi, Eyyi, Eyzi, Ezzi;
double r3inv, r5inv, r7inv, r9inv, rc3inv;
double Ucoul, fcoul, Ucoul_shift;
double smooth_factor;
double U_ind, f_ind[3];

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

int index;
linkedcell->calculateIndex(r_probe,index);

if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
}

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;

int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r = sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz; 

     r3inv = rinv*r2inv;
     r5inv = r3inv*r2inv;
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = qj*r3inv*smooth_factor;
     qr5inv = -3.0*qj*r5inv*smooth_factor;
     qrc = -3.0*qj*rc3inv*r2inv; 

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;

    } 	      
  }
}

// calculation induced dipole interaction on the second neighbors cells
neighborscells = linkedcell->neighbors2_cells[index];
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors2_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r = sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz;

     r3inv = rinv*r2inv;
     r5inv = r3inv*r2inv;
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = qj*r3inv*smooth_factor;
     qr5inv = -3.0*qj*r5inv*smooth_factor;
     qrc = -3.0*qj*rc3inv*r2inv; 

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;
    }
  }
}

f[0] = fx + alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
f[1] = fy + alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
f[2] = fz + alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

Up = U - 0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
return;
}

void Force::coulomb_induced_dipole_aniso_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, U;
double dx, dy, dz;
double r2, x2, y2, z2, r;
double r2inv, r6inv, rinv;
double qi, qj, qr3inv, qr5inv, qrc;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double Exi, Eyi, Ezi, Exxi, Exyi, Exzi, Eyyi, Eyzi, Ezzi;
double r3inv, r5inv, r7inv, r9inv, rc3inv;
double Ucoul, fcoul, Ucoul_shift;
double smooth_factor;
double U_ind, f_ind[3];
double EE[3][3], E[3], alpha_tensor[3][3];
double alpha_parallel, alpha_perpendicular;
double theta, phi;
double un[3];

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];

int index;
linkedcell->calculateIndex(r_probe,index);

if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
}

//polarizability (axial)  = 2.19609742e-30 in m3
//polarizability (radial) = 1.51148405e-30 in m3
alpha_parallel = 2.19609742; 
alpha_perpendicular = 1.51148405;
alpha_parallel *= ALPHA_TO_KCAL_MOL;
alpha_perpendicular *= ALPHA_TO_KCAL_MOL;

// orientation of nitrogen molecule
un[0] = gas->x[0] - gas->x[1];
un[1] = gas->y[0] - gas->y[1];
un[2] = gas->z[0] - gas->z[1];

double un2;
un2 = sqrt(un[0]*un[0] + un[1]*un[1] + un[2]*un[2]);

un[0] /= un2;
un[1] /= un2;
un[2] /= un2;

theta = acos(un[2]);
phi = atan2(un[1],un[0]);

alpha_tensor[0][0] = alpha_perpendicular*pow(sin(phi),2) + pow(cos(phi),2)*(alpha_perpendicular*pow(cos(theta),2)+alpha_parallel*pow(sin(theta),2));
alpha_tensor[1][1] = alpha_perpendicular*pow(cos(phi),2) + pow(sin(phi),2)*(alpha_perpendicular*pow(cos(theta),2)+alpha_parallel*pow(sin(theta),2));
alpha_tensor[2][2] = alpha_parallel*pow(cos(theta),2) + alpha_perpendicular*pow(sin(theta),2);
alpha_tensor[0][1] = (alpha_parallel-alpha_perpendicular)*cos(phi)*sin(phi)*pow(sin(theta),2);
alpha_tensor[1][0] = alpha_tensor[0][1];
alpha_tensor[0][2] = (alpha_parallel-alpha_perpendicular)*cos(phi)*cos(theta)*sin(theta);
alpha_tensor[2][0] = alpha_tensor[0][2];
alpha_tensor[1][2] = (alpha_parallel-alpha_perpendicular)*cos(theta)*sin(phi)*sin(theta);
alpha_tensor[2][1] = alpha_tensor[1][2];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;

int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r = sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz; 

     r3inv = rinv*r2inv;
     r5inv = r3inv*r2inv;
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = qj*r3inv*smooth_factor;
     qr5inv = -3.0*qj*r5inv*smooth_factor;
     qrc = -3.0*qj*rc3inv*r2inv; 

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;

    } 	      
  }
}

// calculation induced dipole interaction on the second neighbors cells
neighborscells = linkedcell->neighbors2_cells[index];
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors2_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r = sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz;

     r3inv = rinv*r2inv;
     r5inv = r3inv*r2inv;
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = qj*r3inv*smooth_factor;
     qr5inv = -3.0*qj*r5inv*smooth_factor;
     qrc = -3.0*qj*rc3inv*r2inv; 

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;
    }
  }
}

// induced-dipole contribution
E[0] = Ex;
E[1] = Ey;
E[2] = Ez;

EE[0][0] = Exx; //Exx
EE[0][1] = Exy; //Exy
EE[0][2] = Exz; //Exz

EE[1][0] = Exy; //Eyx
EE[1][1] = Eyy; //Eyy
EE[1][2] = Eyz; //Eyz

EE[2][0] = Exz; //Ezx
EE[2][1] = Eyz; //Ezy
EE[2][2] = Ezz; //Ezz

U_ind = 0.0;
for (int i = 0; i < 3; i++) {
  for (int j = 0; j < 3; j++) {
    U_ind += -0.5*alpha_tensor[i][j]*E[j]*E[i];
  }
}
Up = U + U_ind;

for (int k = 0; k < 3; k++) {
  f_ind[k] = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      f_ind[k] += 0.5*alpha_tensor[i][j]*(EE[j][k]*E[i] + E[j]*EE[i][k]);
    }
  }
}
f[0] = fx + f_ind[0];
f[1] = fy + f_ind[1];
f[2] = fz + f_ind[2];
return;
}

/*
 * Compute the lennard jones force and potential using linked-cell for CO2
 */
void Force::lennardjones_LC_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy, fz, U, Ulj, flj, Ulj_cut;    
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, r;
double r2inv, r6inv, rc6inv;
double lj1, lj2, lj3, lj4;
double s1, s2;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;
    
int index; 
linkedcell->calculateIndex(r_probe,index); 
  
if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);	    
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
} 	    

int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;  
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];  
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);

    if (r < lj_cutoff) {
      epsilon_target = moleculeTarget->eps[target_id]; 
      sigma_target = moleculeTarget->sig[target_id]; 

      epsilon = sqrt(epsilon_target * epsilon_probe);
      sigma = 0.5*(sigma_target + sigma_probe);
  
      r2inv = 1.0/r2;
      r6inv = r2inv*r2inv*r2inv;
      lj1 = 4.0*epsilon*pow(sigma,6.0);
      lj2 = lj1*pow(sigma,6.0);
      Ulj = r6inv*(lj2*r6inv - lj1);
      rc6inv = 1.0/pow(lj_cutoff,6);    
      Ulj_cut = rc6inv*(lj2*rc6inv - lj1); 
      lj3 = 6.0*lj1;
      lj4 = 12.0*lj2;
      flj = r6inv*(lj4*r6inv - lj3)*r2inv;

      U += Ulj - Ulj_cut;
      fx += flj*dx;
      fy += flj*dy;
      fz += flj*dz;
    }
  }	    
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;
return;
}

/*
 * Compute the lennard jones force and potential for CO2
 */

void Force::lennardjones_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2,r;
double r2inv,r6inv;
double lj1,lj2,lj3,lj4;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];

Ulj = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r = sqrt(r2);

  epsilon_target = moleculeTarget->eps[i]; 
  sigma_target = moleculeTarget->sig[i]; 

  epsilon = sqrt(epsilon_target * epsilon_probe);
  sigma = 0.5*(sigma_target + sigma_probe);
       
  r2inv = 1.0/r2;
  r6inv = r2inv*r2inv*r2inv;
  lj1 = 4.0*epsilon*pow(sigma,6.0);
  lj2 = lj1*pow(sigma,6.0);
  Ulj += r6inv*(lj2*r6inv - lj1);

  lj3 = 6.0*lj1;
  lj4 = 12.0*lj2;
  flj = r6inv*(lj4*r6inv - lj3)*r2inv;
  fx += flj*dx;
  fy += flj*dy;
  fz += flj*dz;
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = Ulj;

return;
}

/*
 * Compute the lennard jones and coulomb interactions using linked-cell CO2
 */
void Force::lennardjones_coulomb_LC_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2,r;
double rinv,r2inv,r6inv;
double lj1,lj2,lj3,lj4;
double qi, qj;
double Ucoul, fcoul, s1, s2;
int index;
double Ulj_cut, rc6inv, Ucoul_shift;


r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];
epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

linkedcell->calculateIndex(r_probe,index);

if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
}

int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);

    if (r < lj_cutoff) {
     epsilon_target = moleculeTarget->eps[target_id]; 
     sigma_target = moleculeTarget->sig[target_id]; 
     epsilon = sqrt(epsilon_target * epsilon_probe);
     sigma = 0.5*(sigma_target + sigma_probe);
     r2inv = 1.0/r2;
     r6inv = r2inv*r2inv*r2inv;
     lj1 = 4.0*epsilon*pow(sigma,6.0);
     lj2 = lj1*pow(sigma,6.0);
     Ulj = r6inv*(lj2*r6inv - lj1);
     rc6inv = 1.0/pow(lj_cutoff,6);
     Ulj_cut = rc6inv*(lj2*rc6inv - lj1);
     lj3 = 6.0*lj1;
     lj4 = 12.0*lj2;
     flj = r6inv*(lj4*r6inv - lj3)*r2inv;		

     U += Ulj - Ulj_cut;
     fx += flj*dx;
     fy += flj*dy;
     fz += flj*dz;
    }

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));
     	
     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz; 
    } 	      
  }
}

// calculation coulomb interaction on the second neighbors cells
neighborscells = linkedcell->neighbors2_cells[index];
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors2_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3); 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));    

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz;
    }
  }
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;

return;
}


/**
 * Compute the lennard jones and coulomb interactions
 **/
void Force::lennardjones_coulomb_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U, Ucoul, fcoul;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, r;
double rinv, r2inv, r6inv;
double lj1, lj2, lj3, lj4;
double qi, qj;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];
epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  epsilon_target = moleculeTarget->eps[i]; 
  sigma_target = moleculeTarget->sig[i]; 

  epsilon = sqrt(epsilon_target * epsilon_probe);
  sigma = 0.5*(sigma_target + sigma_probe);
     
  r2inv = 1.0/r2; 
  r6inv = r2inv*r2inv*r2inv;
  lj1 = 4.0*epsilon*pow(sigma,6.0);
  lj2 = lj1*pow(sigma,6.0);
  Ulj = r6inv*(lj2*r6inv - lj1);
  lj3 = 6.0*lj1;
  lj4 = 12.0*lj2;
  flj = r6inv*(lj4*r6inv - lj3)*r2inv;
  
  U += Ulj;
  fx += flj*dx;
  fy += flj*dy;
  fz += flj*dz;
  
  qj = moleculeTarget->q[i];
  rinv = 1.0/r;
  Ucoul = qi*qj*rinv*KCOUL;
  r2inv = 1.0/r2;
  fcoul = Ucoul*r2inv; 
  U += Ucoul;
  fx += fcoul*dx;
  fy += fcoul*dy;
  fz += fcoul*dz; 
}

f[0] = fx;
f[1] = fy;
f[2] = fz;
Up = U;

return;
}

/**
 * Compute the lennard jones, coulomb and induced dipole interactions using linked-cell CO2 (<- C of CO2)
 **/
void Force::lennardjones_coulomb_induced_dipole_iso_LC_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U, Ulj_cut;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, x2, y2, z2, r;
double r2inv, r6inv, rc6inv;
double lj1, lj2, lj3, lj4;
double q, qr3inv, qr5inv, qrc;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double r3inv, r5inv, r7inv, r9inv;
double rc3inv, smooth_factor;
double rinv;
double Ucoul, Ucoul_shift;
double qi, qj;
double fcoul;


r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];
epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

int index;
linkedcell->calculateIndex(r_probe,index);

if (index >= linkedcell->Ncells || index < 0) {
  printf("outside cell simulation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  Up = 0.0;
  return;
}

Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;

double s1, s2, Exi, Eyi, Ezi, Exxi, Eyyi, Ezzi, Exyi, Exzi, Eyzi;

// calculation lennard-jones and induced dipole interactions on the first neighbors cells
int neighborscells, neighbors, target_id;
neighborscells = linkedcell->neighbors1_cells[index];
int cell_index;
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors1_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);
    r2inv = 1.0/r2;
    
    // lennard-jones interaction
    if (r < lj_cutoff) {
     epsilon_target = moleculeTarget->eps[target_id]; 
     sigma_target = moleculeTarget->sig[target_id]; 

     epsilon = sqrt(epsilon_target * epsilon_probe);
     sigma = 0.5*(sigma_target + sigma_probe);

     r6inv = r2inv*r2inv*r2inv;
     lj1 = 4.0*epsilon*pow(sigma,6.0);
     lj2 = lj1*pow(sigma,6.0);
     Ulj = r6inv*(lj2*r6inv - lj1);
     rc6inv = 1.0/pow(lj_cutoff,6);
     Ulj_cut = rc6inv*(lj2*rc6inv - lj1);

     lj3 = 6.0*lj1;
     lj4 = 12.0*lj2;
     flj = r6inv*(lj4*r6inv - lj3)*r2inv;
       
     U += Ulj - Ulj_cut;
     fx += flj*dx;
     fy += flj*dy;
     fz += flj*dz;
    }

    // coulomb and ion-induced dipole interaction
    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3);
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz;

     r3inv = 1.0/r*r2inv;
     r5inv = r3inv*r2inv;
     
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = qj*r3inv*smooth_factor;
     qr5inv = -3.0*qj*r5inv*smooth_factor;
     qrc = -3.0*qj*rc3inv*r2inv; 

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;
    }
  }
}

// calculation induced dipole interaction on the second neighbors cells
neighborscells = linkedcell->neighbors2_cells[index];
for (int i = 0; i < neighborscells; i++) {
  cell_index = linkedcell->neighbors2_cells_ids[index][i];
  neighbors = linkedcell->atoms_inside_cell[cell_index];
  #pragma omp simd
  for (int j = 0; j < neighbors; j++) {
    target_id = linkedcell->atoms_ids[cell_index][j];
    dx = r_probe[0] - moleculeTarget->x[target_id];
    dy = r_probe[1] - moleculeTarget->y[target_id];
    dz = r_probe[2] - moleculeTarget->z[target_id];
    r2 = dx*dx + dy*dy + dz*dz;
    r =  sqrt(r2);
    r2inv = 1.0/r2;

    // ion-induced dipole interaction
    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
     rinv = 1.0/r;
     Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = -1.5*qi*qj*KCOUL/coul_cutoff + 0.5*qi*qj*KCOUL*r2/pow(coul_cutoff,3);
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv*(1.0 - pow(r/coul_cutoff,3));

     U += Ucoul + Ucoul_shift;
     fx += fcoul*dx;
     fy += fcoul*dy;
     fz += fcoul*dz;
 
     r3inv = rinv*r2inv;
     r5inv = r3inv*r2inv; 
     rc3inv = 1.0/pow(coul_cutoff,3);
     smooth_factor = (1.0 - r*r2*rc3inv);
     qr3inv = qj*r3inv*smooth_factor;
     qr5inv = -3.0*qj*r5inv*smooth_factor;
     qrc = -3.0*qj*rc3inv*r2inv;

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv + dx*dx*qr5inv + dx*dx*qrc;
     Eyyi = qr3inv + dy*dy*qr5inv + dy*dy*qrc;
     Ezzi = qr3inv + dz*dz*qr5inv + dz*dz*qrc;

     Exyi = dx*dy*qr5inv + dx*dy*qrc;
     Exzi = dx*dz*qr5inv + dx*dz*qrc;
     Eyzi = dy*dz*qr5inv + dy*dz*qrc;

     Ex += Exi;
     Ey += Eyi;
     Ez += Ezi;

     Exx += Exxi;
     Eyy += Eyyi;
     Ezz += Ezzi;

     Exy += Exyi;
     Exz += Exzi;
     Eyz += Eyzi;
    }
  }
}

f[0] = fx + alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
f[1] = fy + alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
f[2] = fz + alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

Up = U - 0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
return;
}

/*
 * Compute the lennard jones, coulomb and induced dipole interactions on CO2
 */

void Force::lennardjones_coulomb_induced_dipole_iso_CO2(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy,fz, Ulj, flj, U, Ucoul, fcoul;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, x2, y2, z2, r;
double r2inv, r6inv;
double lj1, lj2, lj3, lj4;
double qi, qj, qr3inv, qr5inv;
double Ex, Ey, Ez, Exx, Exy, Exz, Eyy, Eyz, Ezz;
double rinv, r3inv, r5inv, r7inv, r9inv;

r_probe[0] = gas->x[iatom];
r_probe[1] = gas->y[iatom];
r_probe[2] = gas->z[iatom];
qi = gas->q[iatom];
epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  epsilon_target = moleculeTarget->eps[i];
  sigma_target = moleculeTarget->sig[i];

  epsilon = sqrt(epsilon_target * epsilon_probe);
  sigma = 0.5*(sigma_target + sigma_probe);

  r2inv = 1.0/r2;
  r6inv = r2inv*r2inv*r2inv;
  lj1 = 4.0*epsilon*pow(sigma,6.0);
  lj2 = lj1*pow(sigma,6.0);
  Ulj = r6inv*(lj2*r6inv - lj1);  
  lj3 = 6.0*lj1;
  lj4 = 12.0*lj2;
  flj = r6inv*(lj4*r6inv - lj3)*r2inv;
  U += Ulj;
  fx += flj*dx;
  fy += flj*dy;
  fz += flj*dz;

  qj = moleculeTarget->q[i];
  rinv = 1.0/r;
  Ucoul = qi*qj*rinv*KCOUL;
  r2inv = 1.0/r2;
  fcoul = Ucoul*r2inv;

  U += Ucoul;
  fx += fcoul*dx;
  fy += fcoul*dy;
  fz += fcoul*dz;

  r3inv = rinv*r2inv;
  r5inv = r3inv*r2inv;
  qr3inv = qj*r3inv;
  qr5inv = qj*r5inv;

  Ex += dx * qr3inv;
  Ey += dy * qr3inv;
  Ez += dz * qr3inv;

  Exx += qr3inv - 3.0*dx*dx*qr5inv;
  Eyy += qr3inv - 3.0*dy*dy*qr5inv;
  Ezz += qr3inv - 3.0*dz*dz*qr5inv;

  Exy += - 3.0 * dx * dy * qr5inv;
  Exz += - 3.0 * dx * dz * qr5inv;
  Eyz += - 3.0 * dy * dz * qr5inv;
}

f[0] = fx + alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
f[1] = fy + alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
f[2] = fz + alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

Up = U - 0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
return;
}