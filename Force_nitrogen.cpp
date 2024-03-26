// Samuel Cajahuaringa

#include "headers/Force.h"

Force::Force(MoleculeTarget *moleculeTarget, LinkedCell *linkedcell, double lj_cutoff, double inner_lj_cutoff, 
double alpha, double coul_cutoff, double inner_coul_cutoff) {
this->moleculeTarget = moleculeTarget; 
this->linkedcell = linkedcell;           
this->lj_cutoff = lj_cutoff;                     
this->inner_lj_cutoff =  inner_lj_cutoff;
this->alpha = alpha;
this->coul_cutoff =  coul_cutoff;
this->inner_coul_cutoff =  inner_coul_cutoff;

}

Force::~Force(){	
}

/*
 * Compute the lennard jones force and potential using linked-cell
 */
void Force::lennardjones_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy, fz, U, Ulj, flj;    
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, r;
double r2inv, r6inv;
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
  
if (index >= linkedcell->Ncells && index < 0) {
  printf("wrong calculation: %i\n",index);	    
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
      epsilon = epsilon_target; //sqrt(epsilon_target * epsilon_probe);
      sigma = sigma_target; //0.5*(sigma_target + sigma_probe);
  
      r2inv = 1.0/r2;
      r6inv = r2inv*r2inv*r2inv;
      lj1 = 4.0*epsilon*pow(sigma,6.0);
      lj2 = lj1*pow(sigma,6.0);
      Ulj = r6inv*(lj2*r6inv - lj1);
      lj3 = 6.0*lj1;
      lj4 = 12.0*lj2;
      flj = r6inv*(lj4*r6inv - lj3)*r2inv;

      if (r > inner_lj_cutoff) {
        s1 = Switch(inner_lj_cutoff, lj_cutoff, r);
        s2 = DSwitch(inner_lj_cutoff, lj_cutoff, r);
        flj = s1*flj - s2*Ulj;
        Ulj *= s1;
      }

      U += Ulj;
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

epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

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
  epsilon_target = moleculeTarget->eps[i];
  sigma_target = moleculeTarget->sig[i];
  epsilon = epsilon_target; //sqrt(epsilon_target * epsilon_probe);
  sigma = sigma_target; //0.5*(sigma_target + sigma_probe);
  //epsilon = sqrt(epsilon_target * epsilon_probe);
  //sigma = (sigma_target + sigma_probe)/2.0;
     
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

/*
 * Compute the lennard jones and coulomb interactions using linked-cell
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
double Ulj_shift, rc6inv, Ucoul_shift;


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

if (index >= linkedcell->Ncells && index < 0) {
  printf("wrong calculation: %i\n",index);
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
     //epsilon_target = moleculeTarget->eps[target_id];
     //sigma_target = moleculeTarget->sig[target_id];
     //epsilon = sqrt(epsilon_target * epsilon_probe);
     //sigma = 0.5*(sigma_target + sigma_probe);
     epsilon = moleculeTarget->eps[target_id];
     sigma = moleculeTarget->sig[target_id]; 
     r2inv = 1.0/r2;
     r6inv = r2inv*r2inv*r2inv;
     lj1 = 4.0*epsilon*pow(sigma,6.0);
     lj2 = lj1*pow(sigma,6.0);
     Ulj = r6inv*(lj2*r6inv - lj1);
     rc6inv = 1.0/pow(lj_cutoff,6);
     Ulj_shift = rc6inv*(lj2*rc6inv - lj1);
     lj3 = 6.0*lj1;
     lj4 = 12.0*lj2;
     flj = r6inv*(lj4*r6inv - lj3)*r2inv;

	   /*if (r > inner_lj_cutoff) {
      s1 = Switch(inner_lj_cutoff, lj_cutoff, r);
      s2 = DSwitch(inner_lj_cutoff, lj_cutoff, r);
      flj = s1*flj - s2*Ulj;
      Ulj *= s1; 
     }*/		

     U += Ulj - Ulj_shift;
     fx += flj*dx;
     fy += flj*dy;
     fz += flj*dz;
    }

    if (r < coul_cutoff) {
     qj = moleculeTarget->q[target_id];
	   rinv = 1.0/r;
	   Ucoul = qi*qj*rinv*KCOUL;
     Ucoul_shift = qi*qj*KCOUL/coul_cutoff; 
	   r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv;
     
     /*if (r > inner_coul_cutoff) {
      s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
      s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
      fcoul = s1*fcoul - s2*Ucoul;
      Ucoul *= s1;
     }*/		

     U += Ucoul - Ucoul_shift;
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
     Ucoul_shift = qi*qj*KCOUL/coul_cutoff; 
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv;

     /*if (r > inner_coul_cutoff) {
      s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
      s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
      fcoul = s1*fcoul - s2*Ucoul;
      Ucoul *= s1;
     }*/

     U += Ucoul - Ucoul_shift;
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

epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

double Ucoul, fcoul, s1, s2;

//Ulj = 0.0;
//Ucoul = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  dx = r_probe[0] - moleculeTarget->x[i];
  dy = r_probe[1] - moleculeTarget->y[i];
  dz = r_probe[2] - moleculeTarget->z[i];

  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  //epsilon_target = moleculeTarget->eps[i];
  //sigma_target = moleculeTarget->sig[i];
  //epsilon = sqrt(epsilon_target * epsilon_probe);
  //sigma = 0.5*(sigma_target + sigma_probe);
  //if (r < lj_cutoff) {
  epsilon = moleculeTarget->eps[i];
  sigma = moleculeTarget->sig[i];
     
  r2inv = 1.0/r2; // <-
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
	fcoul = Ucoul*r2inv; // <-

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
 * Compute coulomb interaction using linked-cell
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
  printf("wrong calculation: %i\n",index);
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
     Ucoul_shift = qi*qj*KCOUL/coul_cutoff; 
	   r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv;

	   /*if (r > inner_coul_cutoff) {
      s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
      s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
      fcoul = s1*fcoul - s2*Ucoul;
	    Ucoul *= s1;
     }*/		

     U += Ucoul - Ucoul_shift;
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
     Ucoul_shift = qi*qj*KCOUL/coul_cutoff;      
     r2inv = 1.0/r2;
     fcoul = Ucoul*r2inv;

     /*if (r > inner_coul_cutoff) {
      s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
      s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
      fcoul = s1*fcoul - s2*Ucoul;
      Ucoul *= s1;
     }*/

     U += Ucoul - Ucoul_shift;
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
 * Compute the lennard jones and induced dipole interactions using linked-cell
 **/
void Force::lennardjones_induced_dipole_LC(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
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

epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

U = 0.0;
fx = 0.0;
fy = 0.0;
fz = 0.0;

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

     lj3 = 6.0*lj1;
     lj4 = 12.0*lj2;
     flj = r6inv*(lj4*r6inv - lj3)*r2inv;
       
     if (r > inner_lj_cutoff) {
      s1 = Switch(inner_lj_cutoff, lj_cutoff, r);
      s2 = DSwitch(inner_lj_cutoff, lj_cutoff, r);
      flj = s1*flj - s2*Ulj;
      Ulj *= s1; 
     }

     U += Ulj;
     fx += flj*dx;
     fy += flj*dy;
     fz += flj*dz;
    }

    // ion-induced dipole interaction
    if (r < coul_cutoff) {
     r3inv = 1.0/r*r2inv;
     r5inv = r3inv*r2inv;
     q = moleculeTarget->q[target_id];
     qr3inv = q*r3inv;
     qr5inv = q*r5inv;

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv - 3.0*dx*dx*qr5inv;
     Eyyi = qr3inv - 3.0*dy*dy*qr5inv;
     Ezzi = qr3inv - 3.0*dz*dz*qr5inv;

     Exyi = - 3.0 * dx * dy * qr5inv;
     Exzi = - 3.0 * dx * dz * qr5inv;
     Eyzi = - 3.0 * dy * dz * qr5inv;

     if (r > inner_coul_cutoff) {
      s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
      s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);

      Exxi = s2*dx*Exi + s1*Exxi;
      Eyyi = s2*dy*Eyi + s1*Eyyi;
      Ezzi = s2*dz*Ezi + s1*Ezzi;

      Exyi = s2*dx*Eyi + s1*Exyi;
      Exzi = s2*dx*Ezi + s1*Exzi;
      Eyzi = s2*dy*Ezi + s1*Eyzi;

      Exi *= s1;
      Eyi *= s1;
      Ezi *= s1;
     }

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
     qr3inv = q*r3inv;
     qr5inv = q*r5inv;

     Exi = dx * qr3inv;
     Eyi = dy * qr3inv;
     Ezi = dz * qr3inv;

     Exxi = qr3inv - 3.0*dx*dx*qr5inv;
     Eyyi = qr3inv - 3.0*dy*dy*qr5inv;
     Ezzi = qr3inv - 3.0*dz*dz*qr5inv;

     Exyi = - 3.0 * dx * dy * qr5inv;
     Exzi = - 3.0 * dx * dz * qr5inv;
     Eyzi = - 3.0 * dy * dz * qr5inv;

     if (r > inner_coul_cutoff) {
      s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
      s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);

      Exxi = s2*dx*Exi + s1*Exxi;
      Eyyi = s2*dy*Eyi + s1*Eyyi;
      Ezzi = s2*dz*Ezi + s1*Ezzi;

      Exyi = s2*dx*Eyi + s1*Exyi;
      Exzi = s2*dx*Ezi + s1*Exzi;
      Eyzi = s2*dy*Ezi + s1*Eyzi;

      Exi *= s1;
      Eyi *= s1;
      Ezi *= s1;
     }

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
 * Compute the lennard jones and induced dipole interactions
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

epsilon_probe = gas->eps[iatom];
sigma_probe = gas->sig[iatom];

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

double Force::Switch(double rIn, double rOut, double r) {
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
}

/**
 * Compute the nitrogen interactions using linked-cell
 */
void Force::nitrogen(GasBuffer *gas, vector<double> &f, double &Up) {
double r_probe[3], r_target[3];
double r_N[6][3];
double fx, fy,fz, Ulj, flj, U;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2,r;
double rinv,r2inv,r3inv,r5inv,r6inv;
double lj1,lj2,lj3,lj4;
double qC, qN, q, d;
double Ucoul, fcoul;
int index;

r_probe[0] = gas->x[0];
r_probe[1] = gas->y[0];
r_probe[2] = gas->z[0];
qC = gas->q[0];
qN = -0.5*qC;
d = gas->d;

// Nitrogen 1 axis X
r_N[0][0] = r_probe[0] + gas->x[1];
r_N[0][1] = r_probe[1] + gas->y[1];
r_N[0][2] = r_probe[2] + gas->z[1];
// Nitrogen 2 axis X
r_N[1][0] = r_probe[0] + gas->x[2];
r_N[1][1] = r_probe[1] + gas->y[2];
r_N[1][2] = r_probe[2] + gas->z[2];
// Nitrogen 1 axis Y
r_N[2][0] = r_probe[0] + gas->x[3];
r_N[2][1] = r_probe[1] + gas->y[3];
r_N[2][2] = r_probe[2] + gas->z[3];
// Nitrogen 2 axis Y
r_N[3][0] = r_probe[0] + gas->x[4];
r_N[3][1] = r_probe[1] + gas->y[4];
r_N[3][2] = r_probe[2] + gas->z[4];
// Nitrogen 1 axis Z
r_N[4][0] = r_probe[0] + gas->x[5];
r_N[4][1] = r_probe[1] + gas->y[5];
r_N[4][2] = r_probe[2] + gas->z[5];
// Nitrogen 2 axis Z
r_N[5][0] = r_probe[0] + gas->x[6];
r_N[5][1] = r_probe[1] + gas->y[6];
r_N[5][2] = r_probe[2] + gas->z[6];

double Ulj_N[6], flj_x[6], flj_y[6], flj_z[6];
double Ucoul_N[6], fcoul_Nx[6], fcoul_Ny[6], fcoul_Nz[6];
for (int i = 0; i < 6; i++) {
  Ulj_N[i] = 0.0;
  flj_x[i] = 0.0;
  flj_y[i] = 0.0;
  flj_z[i] = 0.0;
  Ucoul_N[i] = 0.0;
  fcoul_Nx[i] = 0.0;
  fcoul_Ny[i] = 0.0;
  fcoul_Nz[i] = 0.0;
}

double qr3inv, qr5inv;
double Ex, Ey, Ez, Exx, Eyy, Ezz, Exy, Exz, Eyz, Ucoul_C, fcoul_Cx, fcoul_Cy, fcoul_Cz;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;
Ucoul_C = 0.0;
fcoul_Cx = 0.0;
fcoul_Cy = 0.0;
fcoul_Cz = 0.0;

#pragma omp simd
for (int i = 0; i < moleculeTarget->natoms; i++) {
  r_target[0] = moleculeTarget->x[i];
  r_target[1] = moleculeTarget->y[i];
  r_target[2] = moleculeTarget->z[i];

  epsilon = moleculeTarget->eps[i];
  sigma = moleculeTarget->sig[i];
  q = moleculeTarget->q[i];

  // nitrogen calculations
  for (int k = 0; k < 6; k++) {
    dx = r_N[k][0] - r_target[0];
    dy = r_N[k][1] - r_target[1];
    dz = r_N[k][2] - r_target[2];
    r2 = dx*dx + dy*dy + dz*dz;
    r = sqrt(r2);
    
    //if (r < lj_coul) {
 	  r2inv = 1.0/r2;     
    r6inv = r2inv*r2inv*r2inv;
    lj1 = 4.0*epsilon*pow(sigma,6.0);
    lj2 = lj1*pow(sigma,6.0);
    Ulj = r6inv*(lj2*r6inv - lj1);

    lj3 = 6.0*lj1;
    lj4 = 12.0*lj2;
    flj = r6inv*(lj4*r6inv - lj3)*r2inv;

	  Ulj_N[k] += Ulj;
    flj_x[k] += flj*dx;
    flj_y[k] += flj*dy;
    flj_z[k] += flj*dz;
    //}

    //if (r < coul_cutoff) {
	  rinv = 1.0/r;
	  Ucoul = qN*q*rinv*KCOUL;
    r2inv = 1.0/r2;
    fcoul = Ucoul*r2inv;

    Ucoul_N[k] += Ucoul;
	  fcoul_Nx[k] += fcoul*dx;
	  fcoul_Ny[k] += fcoul*dy;
    fcoul_Nz[k] += fcoul*dz; 
    //}
  }

  // central particle calculations
  dx = r_probe[0] - r_target[0];
  dy = r_probe[1] - r_target[1];
  dz = r_probe[2] - r_target[2];
  r2 = dx*dx + dy*dy + dz*dz;
  r =  sqrt(r2);

  //if (r < coul_cutoff) {
	rinv = 1.0/r;
	Ucoul = qC*q*rinv*KCOUL;
	r2inv = 1.0/r2;
  fcoul = Ucoul*r2inv;
  
  Ucoul_C += Ucoul;
	fcoul_Cx += fcoul*dx;
	fcoul_Cy += fcoul*dy;
  fcoul_Cz += fcoul*dz; 
  
  r3inv = rinv*r2inv;
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
  //}
}

double Uind, find[3];
Uind = -0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
find[0] = alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
find[1] = alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
find[2] = alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

double Umol[3], fmolX[3], fmolY[3], fmolZ[3];
Umol[0] = Ulj_N[0] + Ulj_N[1] + Ucoul_N[0] + Ucoul_N[1] + Uind + Ucoul_C;
Umol[1] = Ulj_N[2] + Ulj_N[3] + Ucoul_N[2] + Ucoul_N[3] + Uind + Ucoul_C;
Umol[2] = Ulj_N[4] + Ulj_N[5] + Ucoul_N[4] + Ucoul_N[5] + Uind + Ucoul_C;

/*double consta;
consta = 6.95e-21;

printf("Umol[0]: %g\n",Umol[0]*consta);
printf("Umol[1]: %g\n",Umol[1]*consta);
printf("Umol[2]: %g\n",Umol[2]*consta);*/

double Umin;

Umin = min(Umol[0],min(Umol[1],Umol[2]));

fmolX[0] = flj_x[0] + flj_x[1] + fcoul_Nx[0] + fcoul_Nx[1] + find[0] + fcoul_Cx;
fmolY[0] = flj_y[0] + flj_y[1] + fcoul_Ny[0] + fcoul_Ny[1] + find[1] + fcoul_Cy;
fmolZ[0] = flj_z[0] + flj_z[1] + fcoul_Nz[0] + fcoul_Nz[1] + find[2] + fcoul_Cz;

fmolX[1] = flj_x[2] + flj_x[3] + fcoul_Nx[2] + fcoul_Nx[3] + find[0] + fcoul_Cx;
fmolY[1] = flj_y[2] + flj_y[3] + fcoul_Ny[2] + fcoul_Ny[3] + find[1] + fcoul_Cy;
fmolZ[1] = flj_z[2] + flj_z[3] + fcoul_Nz[2] + fcoul_Nz[3] + find[2] + fcoul_Cy;

fmolX[2] = flj_x[4] + flj_x[5] + fcoul_Nx[4] + fcoul_Nx[5] + find[0] + fcoul_Cx;
fmolY[2] = flj_y[4] + flj_y[5] + fcoul_Ny[4] + fcoul_Ny[5] + find[1] + fcoul_Cy;
fmolZ[2] = flj_z[4] + flj_z[5] + fcoul_Nz[4] + fcoul_Nz[5] + find[2] + fcoul_Cz;

double dU, Z, kBT, T;
dU =0.0;
Z = 0.0;
T = 500.0;
kBT = BOLTZMANN_K * T * J_TO_eV * eV_TO_KCAL_MOL;

for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  Z += exp(-dU/kBT);
}

Up = 0.0;
f[0] = 0.0;
f[1] = 0.0;
f[2] = 0.0;

double w;
/*
for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  w = exp(-dU/kBT)/Z;
  Up += w * Umol[i];
  f[0] += w * fmolX[i];
  f[1] += w * fmolY[i];
  f[2] += w * fmolZ[i];
}*/

int i = 0;
Up = Umol[i];
f[0] = fmolX[i];
f[1] = fmolY[i];
f[2] = fmolZ[i];

//printf("U(r): %g\n",Up*consta);
//printf("fx: %g\n",f[0]*consta*1.0e10);
//printf("fy: %g\n",f[1]*consta*1.0e10);
//printf("fz: %g\n",f[2]*consta*1.0e10);
return;
}

/**
 * Compute the lennard jones force and potential using linked-cell
 */
/*void Force::lennardjones_LC_shift(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
double r_probe[3];
double fx, fy, fz, U, Ulj, flj;    
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, r;
double r2inv, r6inv, s6inv;
double lj1, lj2, lj3, lj4;
double s1, s2, U_shift;

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
  
if (index >= linkedcell->Ncells && index < 0) {
  printf("wrong calculation: %i\n",index);	    
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
      lj3 = 6.0*lj1;
      lj4 = 12.0*lj2;
      flj = r6inv*(lj4*r6inv - lj3)*r2inv;

      s6inv = 1.0/(pow(lj_cutoff,6));
      U_shift = s6inv*(lj2*s6inv - lj1);

      U += Ulj - U_shift;
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
}*/

/**
 * Compute the nitrogen interactions using linked-cell
 */
void Force::nitrogen2_LC(GasBuffer *gas, vector<double> &f, double &Up) {
double r_probe[3], r_target[3];
double r_N[6][3];
double fx, fy, fz, Ulj, flj, U;
double dx, dy, dz;
double epsilon, sigma;
double epsilon_probe, epsilon_target;
double sigma_probe, sigma_target;
double r2, r;
double rinv,r2inv,r3inv,r5inv,r6inv;
double lj1,lj2,lj3,lj4;
double qC, qN, q, d;
double Ucoul, fcoul;
int index;
bool inside;

// central particle
r_probe[0] = gas->x[0];
r_probe[1] = gas->y[0];
r_probe[2] = gas->z[0];
qC = gas->q[0];
qN = -0.5*qC;

linkedcell->calculateIndex(r_probe,index);

// central particle
if (index >= linkedcell->Ncells && index < 0) {
  printf("wrong calculation: %i\n",index);
  printf("position: %g %g %g\n",r_probe[0],r_probe[1],r_probe[2]);
  printf("outside cell\n");
  inside = false;
} else inside = true;

double qr3inv, qr5inv;
double Ex, Ey, Ez, Exx, Eyy, Ezz, Exy, Exz, Eyz, Ucoul_C, fcoul_Cx, fcoul_Cy, fcoul_Cz;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;
Ucoul_C = 0.0;
fcoul_Cx = 0.0;
fcoul_Cy = 0.0;
fcoul_Cz = 0.0;

double s1, s2, Exi, Eyi, Ezi, Exxi, Eyyi, Ezzi, Exyi, Exzi, Eyzi;
int neighborscells, neighbors, target_id;
int cell_index;

if (inside == true) {  
  neighborscells = linkedcell->neighbors1_cells[index];
  for (int i = 0; i < neighborscells; i++) {;
    cell_index = linkedcell->neighbors1_cells_ids[index][i];
    neighbors = linkedcell->atoms_inside_cell[cell_index];
    #pragma omp simd
    for (int j = 0; j < neighbors; j++) {
      target_id = linkedcell->atoms_ids[cell_index][j];
      r_target[0] = moleculeTarget->x[target_id];
      r_target[1] = moleculeTarget->y[target_id];
      r_target[2] = moleculeTarget->z[target_id];
 
      dx = r_probe[0] - r_target[0];
      dy = r_probe[1] - r_target[1];
      dz = r_probe[2] - r_target[2];
      r2 = dx*dx + dy*dy + dz*dz;
      r = sqrt(r2);

      if (r < coul_cutoff) {
        q = moleculeTarget->q[target_id];
	      rinv = 1.0/r;
	      Ucoul = qC*q*rinv*KCOUL;
	      r2inv = 1.0/r2;
        fcoul = Ucoul*r2inv;

        r3inv = rinv*r2inv;
        r5inv = r3inv*r2inv;
        qr3inv = q*r3inv;
        qr5inv = q*r5inv;

        Exi = dx * qr3inv;
        Eyi = dy * qr3inv;
        Ezi = dz * qr3inv;

        Exxi = qr3inv - 3.0*dx*dx*qr5inv;
        Eyyi = qr3inv - 3.0*dy*dy*qr5inv;
        Ezzi = qr3inv - 3.0*dz*dz*qr5inv;

        Exyi = - 3.0 * dx * dy * qr5inv;
        Exzi = - 3.0 * dx * dz * qr5inv;
        Eyzi = - 3.0 * dy * dz * qr5inv;
        
        if (r > inner_coul_cutoff) {
          s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
          s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
          fcoul = s1*fcoul - s2*Ucoul;
          Ucoul *= s1;

          Exxi = s2*dx*Exi + s1*Exxi;
          Eyyi = s2*dy*Eyi + s1*Eyyi;
          Ezzi = s2*dz*Ezi + s1*Ezzi;

          Exyi = s2*dx*Eyi + s1*Exyi;
          Exzi = s2*dx*Ezi + s1*Exzi;
          Eyzi = s2*dy*Ezi + s1*Eyzi;

          Exi *= s1;
          Eyi *= s1;
          Ezi *= s1;
        }

        Ucoul_C += Ucoul;
	      fcoul_Cx += fcoul*dx;
	      fcoul_Cy += fcoul*dy;
        fcoul_Cz += fcoul*dz; 
  
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
  
  neighborscells = linkedcell->neighbors2_cells[index];  
  for (int i = 0; i < neighborscells; i++) {
    cell_index = linkedcell->neighbors2_cells_ids[index][i];
    neighbors = linkedcell->atoms_inside_cell[cell_index];
    #pragma omp simd
    for (int j = 0; j < neighbors; j++) {
      target_id = linkedcell->atoms_ids[cell_index][j];
      r_target[0] = moleculeTarget->x[target_id];
      r_target[1] = moleculeTarget->y[target_id];
      r_target[2] = moleculeTarget->z[target_id];

      dx = r_probe[0] - r_target[0];
      dy = r_probe[1] - r_target[1];
      dz = r_probe[2] - r_target[2];
      r2 = dx*dx + dy*dy + dz*dz;
      r = sqrt(r2);

      if (r < coul_cutoff) {
        q = moleculeTarget->q[target_id];
	      rinv = 1.0/r;
	      Ucoul = qC*q*rinv*KCOUL;
	      r2inv = 1.0/r2;
        fcoul = Ucoul*r2inv;

        r3inv = rinv*r2inv;
        r5inv = r3inv*r2inv;
        qr3inv = q*r3inv;
        qr5inv = q*r5inv;

        Exi = dx * qr3inv;
        Eyi = dy * qr3inv;
        Ezi = dz * qr3inv;

        Exxi = qr3inv - 3.0*dx*dx*qr5inv;
        Eyyi = qr3inv - 3.0*dy*dy*qr5inv;
        Ezzi = qr3inv - 3.0*dz*dz*qr5inv;

        Exyi = - 3.0 * dx * dy * qr5inv;
        Exzi = - 3.0 * dx * dz * qr5inv;
        Eyzi = - 3.0 * dy * dz * qr5inv;

        if (r > inner_coul_cutoff) {
          s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
          s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
          fcoul = s1*fcoul - s2*Ucoul;
          Ucoul *= s1;

          Exxi = s2*dx*Exi + s1*Exxi;
          Eyyi = s2*dy*Eyi + s1*Eyyi;
          Ezzi = s2*dz*Ezi + s1*Ezzi;

          Exyi = s2*dx*Eyi + s1*Exyi;
          Exzi = s2*dx*Ezi + s1*Exzi;
          Eyzi = s2*dy*Ezi + s1*Eyzi;

          Exi *= s1;
          Eyi *= s1;
          Ezi *= s1;
        }

        Ucoul_C += Ucoul;
	      fcoul_Cx += fcoul*dx;
	      fcoul_Cy += fcoul*dy;
        fcoul_Cz += fcoul*dz; 
  
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
}

// nitrogen particles

// Nitrogen 1 axis X
r_N[0][0] = r_probe[0] + gas->x[1];
r_N[0][1] = r_probe[1] + gas->y[1];
r_N[0][2] = r_probe[2] + gas->z[1];
// Nitrogen 2 axis X
r_N[1][0] = r_probe[0] + gas->x[2];
r_N[1][1] = r_probe[1] + gas->y[2];
r_N[1][2] = r_probe[2] + gas->z[2];
// Nitrogen 1 axis Y
r_N[2][0] = r_probe[0] + gas->x[3];
r_N[2][1] = r_probe[1] + gas->y[3];
r_N[2][2] = r_probe[2] + gas->z[3];
// Nitrogen 2 axis Y
r_N[3][0] = r_probe[0] + gas->x[4];
r_N[3][1] = r_probe[1] + gas->y[4];
r_N[3][2] = r_probe[2] + gas->z[4];
// Nitrogen 1 axis Z
r_N[4][0] = r_probe[0] + gas->x[5];
r_N[4][1] = r_probe[1] + gas->y[5];
r_N[4][2] = r_probe[2] + gas->z[5];
// Nitrogen 2 axis Z
r_N[5][0] = r_probe[0] + gas->x[6];
r_N[5][1] = r_probe[1] + gas->y[6];
r_N[5][2] = r_probe[2] + gas->z[6];

double Ulj_N[6], flj_x[6], flj_y[6], flj_z[6];
double Ucoul_N[6], fcoul_Nx[6], fcoul_Ny[6], fcoul_Nz[6];
double pos_N[3];

for (int k = 0; k < 6; k++) {
  pos_N[0] = r_N[k][0];
  pos_N[1] = r_N[k][1];
  pos_N[2] = r_N[k][2];
  
  linkedcell->calculateIndex(pos_N,index);

  // central particle
  if (index >= linkedcell->Ncells && index < 0) {
    printf("wrong calculation: %i\n",index);
    printf("position: %g %g %g\n",pos_N[0],pos_N[1],pos_N[2]);
    printf("outside cell\n");
    inside = false;
  } else inside = true;
  
  Ulj_N[k] = 0.0;
  flj_x[k] = 0.0;
  flj_y[k] = 0.0;
  flj_z[k] = 0.0;
  Ucoul_N[k] = 0.0;
  fcoul_Nx[k] = 0.0;
  fcoul_Ny[k] = 0.0;
  fcoul_Nz[k] = 0.0;

  if (inside == true) {
    neighborscells = linkedcell->neighbors1_cells[index];
    for (int i = 0; i < neighborscells; i++) {
      cell_index = linkedcell->neighbors1_cells_ids[index][i];
      neighbors = linkedcell->atoms_inside_cell[cell_index];
      #pragma omp simd
      for (int j = 0; j < neighbors; j++) {
        target_id = linkedcell->atoms_ids[cell_index][j];
        r_target[0] = moleculeTarget->x[target_id];
        r_target[1] = moleculeTarget->y[target_id];
        r_target[2] = moleculeTarget->z[target_id];
    
        dx = pos_N[0] - r_target[0];
        dy = pos_N[1] - r_target[1];
        dz = pos_N[2] - r_target[2];
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

          lj3 = 6.0*lj1;
          lj4 = 12.0*lj2;
          flj = r6inv*(lj4*r6inv - lj3)*r2inv;

	        if (r > inner_lj_cutoff) {
            s1 = Switch(inner_lj_cutoff, lj_cutoff, r);
            s2 = DSwitch(inner_lj_cutoff, lj_cutoff, r);
            flj = s1*flj - s2*Ulj;
            Ulj *= s1; 
          }		
          Ulj_N[k] += Ulj;
          flj_x[k] += flj*dx;
          flj_y[k] += flj*dy;
          flj_z[k] += flj*dz;
        }  

        if (r < coul_cutoff) {
          q = moleculeTarget->q[target_id];
	        rinv = 1.0/r;
	        Ucoul = qN*q*rinv*KCOUL;
	        r2inv = 1.0/r2;
          fcoul = Ucoul*r2inv;

          if (r > inner_coul_cutoff) {
            s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
            s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
            fcoul = s1*fcoul - s2*Ucoul;
            Ucoul *= s1;
          }
          Ucoul_N[k] += Ucoul;
	        fcoul_Nx[k] += fcoul*dx;
	        fcoul_Ny[k] += fcoul*dy;
          fcoul_Nz[k] += fcoul*dz; 
        }
      }
    }
    
    neighborscells = linkedcell->neighbors2_cells[index];
    for (int i = 0; i < neighborscells; i++) {
      cell_index = linkedcell->neighbors2_cells_ids[index][i];
      neighbors = linkedcell->atoms_inside_cell[cell_index];
      #pragma omp simd
      for (int j = 0; j < neighbors; j++) {
        target_id = linkedcell->atoms_ids[cell_index][j];
        r_target[0] = moleculeTarget->x[target_id];
        r_target[1] = moleculeTarget->y[target_id];
        r_target[2] = moleculeTarget->z[target_id];

        dx = pos_N[0] - r_target[0];
        dy = pos_N[1] - r_target[1];
        dz = pos_N[2] - r_target[2];
        r2 = dx*dx + dy*dy + dz*dz;
        r =  sqrt(r2);  
         
        if (r < coul_cutoff) {
          q = moleculeTarget->q[target_id];
	        rinv = 1.0/r;
	        Ucoul = qN*q*rinv*KCOUL;
	        r2inv = 1.0/r2;
          fcoul = Ucoul*r2inv;
         
          if (r > inner_coul_cutoff) {
            s1 = Switch(inner_coul_cutoff, coul_cutoff, r);
            s2 = DSwitch(inner_coul_cutoff, coul_cutoff, r);
            fcoul = s1*fcoul - s2*Ucoul;
            Ucoul *= s1;
          }		
          Ucoul_N[k] += Ucoul;
	        fcoul_Nx[k] += fcoul*dx;
	        fcoul_Ny[k] += fcoul*dy;
          fcoul_Nz[k] += fcoul*dz; 
        } 
      }  
    }  
  }
}

double Uind, find[3];
Uind = -0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);
find[0] = alpha * (Ex*Exx + Ey*Exy + Ez*Exz);
find[1] = alpha * (Ex*Exy + Ey*Eyy + Ez*Eyz);
find[2] = alpha * (Ex*Exz + Ey*Eyz + Ez*Ezz);

double Umol[3], fmolX[3],fmolY[3],fmolZ[3];
Umol[0] = Ulj_N[0] + Ulj_N[1] + Ucoul_N[0] + Ucoul_N[1] + Uind + Ucoul_C;
Umol[1] = Ulj_N[2] + Ulj_N[3] + Ucoul_N[2] + Ucoul_N[3] + Uind + Ucoul_C;
Umol[2] = Ulj_N[4] + Ulj_N[5] + Ucoul_N[4] + Ucoul_N[5] + Uind + Ucoul_C;

/*double consta;
consta = 6.95e-21;

printf("Umol[0]: %g\n",Umol[0]*consta);
printf("Umol[1]: %g\n",Umol[1]*consta);
printf("Umol[2]: %g\n",Umol[2]*consta);*/

double Umin;

Umin = min(Umol[0],min(Umol[1],Umol[2]));

fmolX[0] = flj_x[0] + flj_x[1] + fcoul_Nx[0] + fcoul_Nx[1] + find[0] + fcoul_Cx;
fmolY[0] = flj_y[0] + flj_y[1] + fcoul_Ny[0] + fcoul_Ny[1] + find[1] + fcoul_Cy;
fmolZ[0] = flj_z[0] + flj_z[1] + fcoul_Nz[0] + fcoul_Nz[1] + find[2] + fcoul_Cz;

fmolX[1] = flj_x[2] + flj_x[3] + fcoul_Nx[2] + fcoul_Nx[3] + find[0] + fcoul_Cx;
fmolY[1] = flj_y[2] + flj_y[3] + fcoul_Ny[2] + fcoul_Ny[3] + find[1] + fcoul_Cy;
fmolZ[1] = flj_z[2] + flj_z[3] + fcoul_Nz[2] + fcoul_Nz[3] + find[2] + fcoul_Cy;

fmolX[2] = flj_x[4] + flj_x[5] + fcoul_Nx[4] + fcoul_Nx[5] + find[0] + fcoul_Cx;
fmolY[2] = flj_y[4] + flj_y[5] + fcoul_Ny[4] + fcoul_Ny[5] + find[1] + fcoul_Cy;
fmolZ[2] = flj_z[4] + flj_z[5] + fcoul_Nz[4] + fcoul_Nz[5] + find[2] + fcoul_Cz;

double dU, Z, kBT, T;
Z = 0.0;
T = 500.0;
kBT = BOLTZMANN_K * T * J_TO_eV * eV_TO_KCAL_MOL;

for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  Z += exp(-dU/kBT);
}

Up = 0.0;
f[0] = 0.0;
f[1] = 0.0;
f[2] = 0.0;

/*double w;
for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  w = exp(-dU/kBT)/Z;
  Up += w * Umol[i];
  f[0] += w * fmolX[i];
  f[1] += w * fmolY[i];
  f[2] += w * fmolZ[i];
}*/

int i = 2;
Up = Umol[i];
f[0] = fmolX[i];
f[1] = fmolY[i];
f[2] = fmolZ[i];

//printf("U(r): %g\n",Up*consta);

return;
}

/*
 * Compute the lennard jones and coulomb interactions using linked-cell
 */
void Force::lennardjones_coulomb_LC_shift(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
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
double Ulj_shift, rc6inv, Ucoul_shift;


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

if (index >= linkedcell->Ncells && index < 0) {
  printf("wrong calculation: %i\n",index);
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
     //epsilon_target = moleculeTarget->eps[target_id];
     //sigma_target = moleculeTarget->sig[target_id];
     //epsilon = sqrt(epsilon_target * epsilon_probe);
     //sigma = 0.5*(sigma_target + sigma_probe);
     epsilon = moleculeTarget->eps[target_id];
     sigma = moleculeTarget->sig[target_id]; 
     r2inv = 1.0/r2;
     r6inv = r2inv*r2inv*r2inv;
     lj1 = 4.0*epsilon*pow(sigma,6.0);
     lj2 = lj1*pow(sigma,6.0);
     Ulj = r6inv*(lj2*r6inv - lj1);
     rc6inv = 1.0/pow(lj_cutoff,6);
     Ulj_shift = rc6inv*(lj2*rc6inv - lj1);
     lj3 = 6.0*lj1;
     lj4 = 12.0*lj2;
     flj = r6inv*(lj4*r6inv - lj3)*r2inv;		

     U += Ulj - Ulj_shift;
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
 * Compute coulomb interaction using linked-cell
 **/
void Force::coulomb_LC_shift(GasBuffer *gas, int iatom, vector<double> &f, double &Up) {
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
  printf("wrong calculation: %i\n",index);
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
