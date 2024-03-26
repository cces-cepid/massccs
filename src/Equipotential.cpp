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

#include "headers/Equipotential.h"

Equipotential::Equipotential(MoleculeTarget *moleculeTarget, GasBuffer *gas, unsigned int polarizability_flag, double temperature, double mu, double alpha, unsigned int gas_buffer_flag) {
this->moleculeTarget = moleculeTarget;
this->gas = gas;
this->polarizability_flag = polarizability_flag;
this->temperature = temperature;
this->mu = mu;
this->alpha = alpha;
this->gas_buffer_flag = gas_buffer_flag;

a = 0.0;
b = 0.0;
c = 0.0;

v_mb = 2.0*sqrt((2.0 * BOLTZMANN_K * temperature)/(M_PI * mu * AMU_TO_KG)); // v(m/s)
v2_mb = 3.0 * BOLTZMANN_K * temperature/(mu * AMU_TO_KG); 
sig_mb = sqrt(v2_mb - v_mb*v_mb);
v_min_mb = v_mb - 2.0*sig_mb;
Ek_min = 0.5 * mu * AMU_TO_KG * pow(v_min_mb,2.0) * J_TO_eV * eV_TO_KCAL_MOL;  

maxX = moleculeTarget->maxX;
maxY = moleculeTarget->maxY;
maxZ = moleculeTarget->maxZ;

equipotentialpoints();

ellipsoid();

print();
}

Equipotential::~Equipotential() {
  cout << "equipotential destroyed" << endl;
}

/* this routine is adapted from the CoSIMS software
 * for more information review the following github project
 * https://github.com/ChristopherAMyers/CoSIMS
 */  

void Equipotential::equipotentialpoints() {
/*
This functions generates the potential energy boundry
that will be fit to an ellipsoid. Any particle outside this
boundry (and thus the ellipsoid) is defined to have a
potential energy deemed to be zero (by approximation).
The zero potential energy is chosen such that the kinetic
energy of the slowest possible moving gas particle is only
lowered by 1%. The boundry thus defines the starting and
ending position of all gas particle trajectories.
*/

vector<vector<double>> boundryPoints_temp{};
int numRotations = 20, count = 0;                   //number of rotations about z-axis
double theta;                                       //holds angle of rotation
vector<double> vel(3), pos(3);
vector<double> prevPos(3), blank{0.0,0.0,0.0}; //vectors for determining test points
double velMag = 0.1;                                //incrimental vector
double posMag = abs(max(abs(maxX), abs(maxY)));     //starting distance from z-axis
double startX, startY;                              //start positions
double zPos = 0.0, zMax = 0.0;                      //distance along z-axis, maximum extent along z-axis
double energy;                                      //holds potential at pos
double slope = 1.0;                                 //determins search directions
double zInc = 0.0; int maxZLevels = 25;             //z-axis incrimental values
bool accept = false, hit = false;                   //allows exit from while-loops

//first, find maximum extent along z-axis
accept = false;
pos[0] = 0.0;
pos[1] = 0.0;
pos[2] = abs(maxZ) + 20.0; // initial position equal maxZ + 20 Angstrom
//energy = potential(pos);

if (gas_buffer_flag == 1 || gas_buffer_flag == 4 || gas_buffer_flag == 5) {
  energy = potential_He(pos);
} else if (gas_buffer_flag == 2) {
  energy = potential_N2(pos);
} else if (gas_buffer_flag == 3) {
  energy = potential_CO2(pos);
}

Ek_min = 0.01*Ek_min;

if(energy + Ek_min < 0.0) slope = -1.0; // total energy is negative

while (!accept) {
  prevPos = pos;
  pos[2] -= 0.5*slope; // increment or decrement 0.5 Ansgtrom
  //energy = potential(pos);
  if (gas_buffer_flag == 1 || gas_buffer_flag == 4 || gas_buffer_flag == 5) {
    energy = potential_He(pos);
  } else if (gas_buffer_flag == 2) {
    energy = potential_N2(pos);
  } else if (gas_buffer_flag == 3) {
    energy = potential_CO2(pos);
  }

  //found closest distance along z-axis
  if ((energy + Ek_min)*slope < 0.0) {
    zMax = max(pos[2], prevPos[2]);
    accept = true;
  }
}

//start search from bottom-up
zMax *=1.3;
posMag += abs(zMax - abs(maxZ));
maxZLevels = ceil(max(25, (int)ceil(zMax)));
zPos = -zMax;
zInc = 2*zMax/(maxZLevels);

boundryPoints_temp.resize((int)(maxZLevels + 1)*numRotations, blank);
#pragma omp simd
for (int i = 0; i < (maxZLevels + 1); i ++) {
  zPos = -zMax + i*zInc;
  //loop through each rotation about z-axis
  for (int n = 0; n < numRotations; n ++) {
    //create starting positions
    slope = 1.0;
    theta = n*2.0*M_PI/((double)numRotations);
    startX = posMag*cos(theta);
    startY = posMag*sin(theta);
    pos[0] = startX;
    pos[1] = startY;
    pos[2] = zPos;

    vel[0] = -velMag*cos(theta);
    vel[1] = -velMag*sin(theta);
    vel[2] = 0.0;
    hit = false;

    //get potential energy at pos
    //energy = potential(pos);
    if (gas_buffer_flag == 1 || gas_buffer_flag == 4 || gas_buffer_flag == 5) {
      energy = potential_He(pos);
    } else if (gas_buffer_flag == 2) {
      energy = potential_N2(pos);
    } else if (gas_buffer_flag == 3) {
      energy = potential_CO2(pos);
    }

    //if within the boundry that is being set up,
    //search from inside-out instead
    if(energy + Ek_min < 0.0) slope = -1.0;
    vel[0] = vel[0]*slope;
    vel[1] = vel[1]*slope;
    vel[2] = vel[2]*slope;

    while (!hit) {
      //energy = potential(pos);
      if (gas_buffer_flag == 1 || gas_buffer_flag == 4 || gas_buffer_flag == 5) {
        energy = potential_He(pos);
      } else if (gas_buffer_flag == 2) {
        energy = potential_N2(pos);
      } else if (gas_buffer_flag == 3) {
        energy = potential_CO2(pos);
      }

      //potential energy vs kinetic energy condition
      if ((energy + Ek_min)*slope <  0.0) {
        hit = true;
        //choose point with farthest distance
        if (Math::vecModulus(prevPos) > Math::vecModulus(pos)) {
          boundryPoints_temp[i*numRotations + n][0] = prevPos[0];
          boundryPoints_temp[i*numRotations + n][1] = prevPos[1];
          boundryPoints_temp[i*numRotations + n][2] = prevPos[2];
        } else {
          boundryPoints_temp[i*numRotations + n][0] = pos[0];
          boundryPoints_temp[i*numRotations + n][1] = pos[1];
          boundryPoints_temp[i*numRotations + n][2] = pos[2];
        }
        count++;
      //if energy condition was not met,
      //but traversed across other side of molecule  
      } else if ((startX == 0 && startY/pos[1] <= 0)
        || (startY == 0 && startX/pos[0] <= 0)
        || (startX/pos[0] <= 0 && startY/pos[1] <= 0)) {
          hit = true;
      //update position vectors        
      } else {
        prevPos = pos;
        pos = Math::addVec(pos,vel);
      }
    }
  }
}

boundryPoints.reserve(count);
count = 0;
bool test;
for (int i = 0; i < (signed)boundryPoints_temp.size(); i ++) {
   test = !(boundryPoints_temp[i][0] == 0 && boundryPoints_temp[i][1] == 0 && boundryPoints_temp[i][2] == 0);
   if (!(boundryPoints_temp[i][0] == 0 && boundryPoints_temp[i][1] == 0 && boundryPoints_temp[i][2] == 0)) {
     boundryPoints.push_back(boundryPoints_temp[i]);
     count++;
   }
}

}

double Equipotential::potential_He(vector<double> pos) {
double pot;
double epsilon_target, sigma_target;
double Ulj, Uind;
double Ex, Ey, Ez;

Ulj = 0.0;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;

#pragma omp parallel for reduction(+:Ulj,Ex,Ey,Ez) 
for (unsigned int i = 0; i < moleculeTarget->natoms; i++) {
  double dx = pos[0] - moleculeTarget->x[i];
  double dy = pos[1] - moleculeTarget->y[i];
  double dz = pos[2] - moleculeTarget->z[i];

  double r2 = dx*dx + dy*dy + dz*dz;
  double r2inv = 1.0/r2;
  double r6inv = r2inv*r2inv*r2inv;

  double epsilon = moleculeTarget->eps[i];
  double sigma = moleculeTarget->sig[i];
 
  double lj1 = 4.0*epsilon*pow(sigma,6);
  double lj2 = lj1*pow(sigma,6);
  Ulj += r6inv*(lj2*r6inv - lj1);

  // ion-induced dipole interaction
  if (polarizability_flag != 0) {
    double rinv  = sqrt(r2inv);
    double r3inv = rinv*r2inv;
    double q = moleculeTarget->q[i];
    Ex += q * dx * r3inv;
    Ey += q * dy * r3inv;
    Ez += q * dz * r3inv;
  }  
}

pot = Ulj;

if (polarizability_flag != 0) pot -= 0.5*alpha*(Ex*Ex + Ey*Ey + Ez*Ez);

return pot;
}

double Equipotential::potential_N2(vector<double> pos) {
double pot;
double epsilon_target, sigma_target;
double r_N[6][3], Ulj_N[6], Ucoul_N[6], Ucoul_C;
double Ex, Ey, Ez;
double Exx, Eyy, Ezz, Exy, Exz, Eyz;
double qC, qN, d;

qC = 0.965;
qN = -0.5*qC;
d = 0.5488;

// Nitrogen 1 axis X
r_N[0][0] = pos[0] + d;
r_N[0][1] = pos[1];
r_N[0][2] = pos[2];
// Nitrogen 2 axis X
r_N[1][0] = pos[0] - d;
r_N[1][1] = pos[1];
r_N[1][2] = pos[2];
// Nitrogen 1 axis Y
r_N[2][0] = pos[0];
r_N[2][1] = pos[1] + d;
r_N[2][2] = pos[2];
// Nitrogen 2 axis Y
r_N[3][0] = pos[0];
r_N[3][1] = pos[1] - d;
r_N[3][2] = pos[2];
// Nitrogen 1 axis Z
r_N[4][0] = pos[0];
r_N[4][1] = pos[1];
r_N[4][2] = pos[2] + d;
// Nitrogen 2 axis Z
r_N[5][0] = pos[0];
r_N[5][1] = pos[1];
r_N[5][2] = pos[2] - d;

for (int i = 0; i < 6; i++) {
  Ulj_N[i] = 0.0;
  Ucoul_N[i] = 0.0;
}

Ucoul_C = 0.0;
Ex = 0.0;
Ey = 0.0;
Ez = 0.0;
Exx = 0.0;
Eyy = 0.0;
Ezz = 0.0;
Exy = 0.0;
Exz = 0.0;
Eyz = 0.0;

#pragma omp parallel for reduction(+:Ulj_N,Ucoul_N,Ucoul_C,Ex,Ey,Ez) 
for (int i = 0; i < moleculeTarget->natoms; i++) {
  double r_target[3];
  r_target[0] = moleculeTarget->x[i];
  r_target[1] = moleculeTarget->y[i];
  r_target[2] = moleculeTarget->z[i];
  double q = moleculeTarget->q[i];
  double epsilon = moleculeTarget->eps[i];
  double sigma = moleculeTarget->sig[i];
  double lj1 = 4.0*epsilon*pow(sigma,6);
  double lj2 = lj1*pow(sigma,6);

  // nitrogen calculations
  for (int k =0; k < 6; k++) {
    double dx = r_N[k][0] - r_target[0];
    double dy = r_N[k][1] - r_target[1];
    double dz = r_N[k][2] - r_target[2];

    double r2 = dx*dx + dy*dy + dz*dz;
    double r2inv = 1.0/r2;
    double r6inv = r2inv*r2inv*r2inv;

    double Ulj = r6inv*(lj2*r6inv - lj1);
    Ulj_N[k] += Ulj;

    if (polarizability_flag != 0) {
      double rinv  = sqrt(r2inv);
      double Ucoul = qN*q*rinv*KCOUL;
      Ucoul_N[k] += Ucoul;
    }
  }

  // central particle calculations
  if (polarizability_flag != 0) {
    double dx = pos[0] - r_target[0];
    double dy = pos[1] - r_target[1];
    double dz = pos[2] - r_target[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r2);
    double rinv = 1.0/r;
    double Ucoul = qC*q*rinv*KCOUL;

    Ucoul_C += Ucoul;

    double r2inv = 1.0/r2;
    double r3inv = rinv*r2inv;

    Ex += dx * q * r3inv;
    Ey += dy * q * r3inv;
    Ez += dz * q * r3inv;
  }
}

double Uind;
Uind = -0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);

double Umol[3];
Umol[0] = Ulj_N[0] + Ulj_N[1] + Ucoul_N[0] + Ucoul_N[1] + Uind + Ucoul_C;
Umol[1] = Ulj_N[2] + Ulj_N[3] + Ucoul_N[2] + Ucoul_N[3] + Uind + Ucoul_C;
Umol[2] = Ulj_N[4] + Ulj_N[5] + Ucoul_N[4] + Ucoul_N[5] + Uind + Ucoul_C;

double Umin;
Umin = min(Umol[0],min(Umol[1],Umol[2]));

double dU, Z, kBT, T;
dU =0.0;
Z = 0.0;
T = 500.0;
kBT = BOLTZMANN_K * T * J_TO_eV * eV_TO_KCAL_MOL;

for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  Z += exp(-dU/kBT);
}

double w;
pot = 0.0;
for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  w = exp(-dU/kBT)/Z;
  pot += w * Umol[i];
}

return pot;
}

/* developed in progress */
double Equipotential::potential_CO2(vector<double> pos) {
double pot;
double epsilon_target, sigma_target;
double r_O[6][3], Ulj_O[6], Ucoul_O[6], Ulj_C, Ucoul_C;
double Ex, Ey, Ez;
double Exx, Eyy, Ezz, Exy, Exz, Eyz;
double qC, qO, d;
double eps_O, eps_C;
double sig_O, sig_C;

qC = 0.6512;
qO = -0.5*qC;
d = 1.149;

eps_O = 0.159; 
eps_C = 0.055;   
sig_O = 3.033; 
sig_C = 2.757; 

// Oxygen 1 axis X
r_O[0][0] = pos[0] + d;
r_O[0][1] = pos[1];
r_O[0][2] = pos[2];
// Oxygen 2 axis X
r_O[1][0] = pos[0] - d;
r_O[1][1] = pos[1];
r_O[1][2] = pos[2];
// Oxygen 1 axis Y
r_O[2][0] = pos[0];
r_O[2][1] = pos[1] + d;
r_O[2][2] = pos[2];
// Oxygen 2 axis Y
r_O[3][0] = pos[0];
r_O[3][1] = pos[1] - d;
r_O[3][2] = pos[2];
// Oxygen 1 axis Z
r_O[4][0] = pos[0];
r_O[4][1] = pos[1];
r_O[4][2] = pos[2] + d;
// Oxygen 2 axis Z
r_O[5][0] = pos[0];
r_O[5][1] = pos[1];
r_O[5][2] = pos[2] - d;

Ulj_C = 0.0;
Ucoul_C = 0.0;
for (int i = 0; i < 6; i++) {
  Ulj_O[i] = 0.0;
  Ucoul_O[i] = 0.0;
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

#pragma omp parallel for reduction(+:Ulj_O,Ucoul_O,Ulj_C,Ucoul_C,Ex,Ey,Ez) 
for (int i = 0; i < moleculeTarget->natoms; i++) {
  double r_target[3];
  r_target[0] = moleculeTarget->x[i];
  r_target[1] = moleculeTarget->y[i];
  r_target[2] = moleculeTarget->z[i];
  double q = moleculeTarget->q[i];
  double epsilon = moleculeTarget->eps[i];
  double sigma = moleculeTarget->sig[i];
  double lj1 = 4.0*epsilon*pow(sigma,6);
  double lj2 = lj1*pow(sigma,6);

  // oxygen calculations
  for (int k =0; k < 6; k++) {
    double dx = r_O[k][0] - r_target[0];
    double dy = r_O[k][1] - r_target[1];
    double dz = r_O[k][2] - r_target[2];

    double r2 = dx*dx + dy*dy + dz*dz;
    double r2inv = 1.0/r2;
    double r6inv = r2inv*r2inv*r2inv;
  
    double Ulj = r6inv*(lj2*r6inv - lj1);
    Ulj_O[k] += Ulj;

    if (polarizability_flag != 0) {
      double rinv  = sqrt(r2inv);
      double Ucoul = qO*q*rinv*KCOUL;
      Ucoul_O[k] += Ucoul;
    }
  }

  // carbon particle calculations
  double dx = pos[0] - r_target[0];
  double dy = pos[1] - r_target[1];
  double dz = pos[2] - r_target[2];
  double r2 = dx*dx + dy*dy + dz*dz;
  double r2inv = 1.0/r2;
  double r6inv = r2inv*r2inv*r2inv;

  double epsilon_central = moleculeTarget->eps_central[i];
  double sigma_central = moleculeTarget->sig_central[i];
  double lj1_central= 4.0*epsilon_central*pow(sigma_central,6);
  double lj2_central = lj1_central*pow(sigma_central,6);

  double Ulj = r6inv*(lj2_central*r6inv - lj1_central);
  Ulj_C += Ulj;

  if (polarizability_flag != 0) {    
    double r = sqrt(r2);
    double rinv = 1.0/r;     
    double Ucoul = qC*q*rinv*KCOUL;
    Ucoul_C += Ucoul;

    double r2inv = 1.0/r2;
    double r3inv = rinv*r2inv;

    Ex += dx * q * r3inv;
    Ey += dy * q * r3inv;
    Ez += dz * q * r3inv;
  }
}

double Uind;
Uind = -0.5 * alpha * (Ex*Ex + Ey*Ey + Ez*Ez);

double Umol[3];
Umol[0] = Ulj_O[0] + Ulj_O[1] + Ucoul_O[0] + Ucoul_O[1] + Uind + Ulj_C + Ucoul_C;
Umol[1] = Ulj_O[2] + Ulj_O[3] + Ucoul_O[2] + Ucoul_O[3] + Uind + Ulj_C + Ucoul_C;
Umol[2] = Ulj_O[4] + Ulj_O[5] + Ucoul_O[4] + Ucoul_O[5] + Uind + Ulj_C + Ucoul_C;

double Umin;
Umin = min(Umol[0],min(Umol[1],Umol[2]));

double dU, Z, kBT, T;
dU =0.0;
Z = 0.0;
T = 500.0;
kBT = BOLTZMANN_K * T * J_TO_eV * eV_TO_KCAL_MOL;

for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  Z += exp(-dU/kBT);
}

double w;
pot = 0.0;
for (int i = 0; i < 3; i++) {
  dU = Umol[i] - Umin;
  w = exp(-dU/kBT)/Z;
  pot += w * Umol[i];
}

return pot;
}

void Equipotential::ellipsoid() {
//stores summation information
double mu = 0, nu = 0, epsilon = 0, delta = 0, sigma = 0, rho = 0;
double alpha = 0, beta = 0, gamma = 0;
double A = 0, B = 0, C = 0;
rmsd = 0;

//variables for ellipsoid surface area
double p, f, deriv_a, deriv_b, deriv_c;

//stores uncertainty information
aUncert = 0, bUncert = 0, cUncert = 0;

//double x2 = 0, y2 = 0, z2 = 0, x4 = 0, y4 = 0, z4 = 0, x2y2 = 0, x2z2 = 0, y2z2 = 0;

//calculate sumamtion over data coordinates
#pragma omp parallel for reduction(+:mu,nu,epsilon,delta,sigma,rho,alpha,beta,gamma)
for (int i = 0; i < boundryPoints.size(); i++) {
  double x2 = boundryPoints[i][0]*boundryPoints[i][0];
  double y2 = boundryPoints[i][1]*boundryPoints[i][1];
  double z2 = boundryPoints[i][2]*boundryPoints[i][2];

  double x4 = x2*x2;
  double y4 = y2*y2; 
  double z4 = z2*z2;
  double x2y2 = x2*y2;
  double x2z2 = x2*z2;
  double y2z2 = y2*z2;

  mu = mu + x4;
  nu = nu + y4;
  epsilon = epsilon + z4;

  delta = delta + x2*y2;
  sigma = sigma + x2*z2;
  rho = rho + y2*z2;

  alpha = alpha + x2;
  beta = beta + y2;
  gamma = gamma + z2;
}

double det = epsilon*mu*nu - epsilon*delta*delta - mu*rho*rho + 2 * delta*rho*sigma - nu*sigma*sigma;

A = (alpha*(nu*epsilon - rho*rho)
   + delta*(rho*gamma - beta*epsilon)
   + sigma*(beta*rho - nu*gamma)) / det;

B = (mu*(beta*epsilon - gamma*rho)
   + alpha*(rho*sigma - delta*epsilon)
   + sigma*(delta*gamma - beta*sigma)) / det;

C = (mu*(nu*gamma - beta*rho)
   + delta*(beta*sigma - delta*gamma)
   + alpha*(delta*rho - nu*sigma)) / det;

a = 1.0/sqrt(A); b = 1.0/sqrt(B); c = 1.0/sqrt(C);

rmsd = sqrt((A*A*mu + B*B*nu + C*C*epsilon +
     2 * (A*B*delta + B*C*rho + A*C*sigma - A*alpha - B*beta - C*gamma) + boundryPoints.size()) / (boundryPoints.size() - 1));
cov_aa = 2 * rmsd*rmsd* (epsilon*nu - rho*rho)*a*a     / (4 * A*A * det);
cov_bb = 2 * rmsd*rmsd* (epsilon*mu - sigma*sigma)*b*b / (4 * B*B * det);
cov_cc = 2 * rmsd*rmsd* (mu*nu - delta*delta)*c*c      / (4 * C*C * det);

cov_ab = 2 * rmsd*rmsd* (rho*sigma - delta*epsilon)*a*b / (4 * A*B * det);
cov_ac = 2 * rmsd*rmsd* (delta*rho - nu*sigma)*a*c      / (4 * A*C * det);
cov_bc = 2 * rmsd*rmsd* (delta*sigma - mu*rho)*b*c      / (4 * B*C * det);

aUncert = (a / A)*sqrt((epsilon*nu - rho*rho) * 2 / det)*rmsd/2;
bUncert = (b / B)*sqrt((epsilon*mu - sigma*sigma) * 2 / det)*rmsd/2;
cUncert = (c / C)*sqrt((mu*nu - delta*delta) * 2 / det)*rmsd/2;

avgAxes = (a + b + c) / 3.0;
avgAxesUncert = sqrt(cov_aa + cov_bb + cov_cc + 2 * (cov_ab + cov_ac + cov_bc)) / 3.0;

// aproximate formule of area
p = 1.6075;
f = (pow(a*b, p) + pow(a*c, p) + pow(b*c, p)) / 3;
surfArea = 4 * M_PI*pow(f, 1.0/p);
deriv_a = surfArea*pow(a, p - 1)*(pow(b, p) + pow(c, p)) / (3 * f);
deriv_b = surfArea*pow(b, p - 1)*(pow(a, p) + pow(c, p)) / (3 * f);
deriv_c = surfArea*pow(c, p - 1)*(pow(a, p) + pow(b, p)) / (3 * f);

surfAreaUncert = sqrt(deriv_a*deriv_a*cov_aa + deriv_b*deriv_b*cov_bb + deriv_c*deriv_c*cov_cc
      + 2 * (deriv_a*deriv_b*cov_ab + deriv_a*deriv_c*cov_ac + deriv_b*deriv_c*cov_bc));
}

void Equipotential::print() {

cout << "*********************************************************"
            << endl;
cout << "Equipotential Ellipsoid: " << endl;
cout << "*********************************************************"
            << endl;
cout << "Axis Length: " << a << "  "<< b << "  "<< c << "  Ang" << endl;
cout << "Uncertainy: " << aUncert << "  "<< bUncert << "  "<< cUncert << endl;
cout << "Percent Uncert: " << aUncert*100/a << "  "<< bUncert*100/b << "  "<< cUncert*100/b << endl;
cout << "RMSD: " << rmsd << endl;
cout << "Covariance matrix: " << rmsd << endl;
cout << "{  " << cov_aa << "  " << cov_ab << "  " << cov_ac << "  }" << endl;
cout << "{  " << cov_ab << "  " << cov_bb << "  " << cov_bc << "  }" << endl;
cout << "{  " << cov_ac << "  " << cov_bc << "  " << cov_cc << "  }" << endl;
cout << "Average ellipsoid axes: " << avgAxes << " +/- "<< avgAxesUncert << " Ang" << endl;
cout << "Ellipsoid surface area: " << surfArea << " +/- "<< surfAreaUncert << " Ang^2" << endl;

enlargeEllipsoidBoundry();

cout << "Ellipsoid axes further expanded by: " << enlargeAmount << " Ang" << endl;

}

void Equipotential::enlargeEllipsoidBoundry() {

double aOld = a, bOld = b, cOld = c;

bool accept = false;
double aInc = 0.01;
double bInc = 0.01;
double cInc = 0.01;
int i = 0;

double mag = 0;

int counter = 0;

// continue to enlarge ellipsoid until all
//boundry points are inside
do {
  counter++;
  //loop through all boundry points
  for (i = 0; i < (signed)boundryPoints.size() + 1; i++) {
     //incase running in an infinite loop
     if (i == (signed)boundryPoints.size() || a > 3*aOld) {
       accept = true;
       break;
     }

     //definition of al ellipsoid: (x/z)^2 + (y/b)^2 + (z/c)^2 = 1
     mag = boundryPoints[i][0]*boundryPoints[i][0] / (a*a)
         + boundryPoints[i][1]*boundryPoints[i][1] / (b*b)
         + boundryPoints[i][2]*boundryPoints[i][2] / (c*c);
     if (mag > 1) {
       a = a + aInc;
       b = b + bInc;
       c = c + cInc;
       break;
     }
  }
} while (!accept);

enlargeAmount = a - aOld;

}


