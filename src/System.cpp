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

#include "headers/System.h"

System::System(char *inputFilename) {
// read the input
double start = omp_get_wtime(); 
input = new Input(inputFilename); 

// TODO: move to Output
input->printReadInput(); 

// initialize variables from input values
nProbe = input->nProbe;                           // numbers of gas buffers 
nIter = input->nIter;                             // numbers of CCS calculations
seed = input->seed;                               // random number seed 
mt = new RandomNumber(seed);
nthreads = input->nthreads;                       // number of threads
targetFilename = input->targetFilename;           // xyz or pqr file of molecule target 
dt = input->dt;                                   // time step in fs 
temperatureTarget = input->temperatureTarget;     // temperature in Kelvin
gas_buffer_flag = input->gas_buffer_flag;         // He = 1, N2 = 2
equipotential_flag = input->equipotential_flag; 
skin = input->skin;                               // skin cell-size   
short_range_cutoff = input->short_range_cutoff;   // yes = 1 and not = 0 for cut lennard-jones interacion
lj_cutoff = input->lj_cutoff;                     // lennard-jones cutoff   
long_range_flag =  input->long_range_flag;        // apply coulomb interaction 
long_range_cutoff = input->long_range_cutoff;     // yes = 1 and not = 0 for cut coulomb interacion
coul_cutoff = input->coul_cutoff;                 // coulomb cutoff
polarizability_flag = input->polarizability_flag; // apply induced-dipole interaction
alpha = input->alpha;                             // polarizability
user_ff_flag = input->user_ff_flag;               // user force field
user_ff = input->user_ff;                         // name of force field 

if (gas_buffer_flag == 1 || gas_buffer_flag == 4) {
  if (short_range_cutoff == 0 && long_range_flag == 0) { 
  // only lennard-jones interaction 
  force_type = 1;  
  } else if (short_range_cutoff == 1 && long_range_flag == 0) { 
    // only lennard-jones interaction with cutoff
    force_type = 2;  
  } else if (short_range_cutoff == 0 && polarizability_flag == 1 && long_range_cutoff == 0) {
    // lennard-jones and induced-dipole interactions
    force_type = 3;  
  } else if (short_range_cutoff == 1 && polarizability_flag == 1 && long_range_cutoff == 1) {
    // lennard-jones and induced-dipole interactions with cutoff
    force_type = 4;  
  } else {
  force_type = 2;
  }
} else if (gas_buffer_flag == 2 || gas_buffer_flag == 3){
  if (short_range_cutoff == 0 && long_range_flag == 0 && polarizability_flag == 0) { 
  // only lennard-jones interaction 
  force_type = 1;  
  } else if (short_range_cutoff == 1 && long_range_flag == 0 && polarizability_flag == 0) { 
    // only lennard-jones interaction with cutoff
    force_type = 2;  
  } else if (short_range_cutoff == 0 && long_range_flag == 1 && long_range_cutoff == 0 && polarizability_flag == 0) {
    // lennard-jones and coulomb interactions
    force_type = 3;  
  } else if (short_range_cutoff == 1 && long_range_flag == 1 && long_range_cutoff == 1 && polarizability_flag == 0) {
    // lennard-jones and coulomb interactions with cutoff
    force_type = 4;  
  } else if (short_range_cutoff == 0 && long_range_flag == 1 && long_range_cutoff == 0 && polarizability_flag == 1) {
    // lennard-jones, coulomb and induced-dipole interactions 
    force_type = 5;  
  } else if (short_range_cutoff == 1 && long_range_flag == 1 && long_range_cutoff == 1 && polarizability_flag == 1) {
    // lennard-jones, coulomb and induced-dipole interactions with cutoff
    force_type = 6;  
  } else {
  force_type = 2;
  }
}

gas = new GasBuffer(gas_buffer_flag);

double start_molecule = omp_get_wtime();
moleculeTarget = new MoleculeTarget(targetFilename, gas_buffer_flag, user_ff, user_ff_flag, force_type); 
double end_molecule = omp_get_wtime();
cout << "orientation time of molecule target: " << (end_molecule - start_molecule) << " s" << endl;

// reduced mass
if (gas_buffer_flag == 1 || gas_buffer_flag == 4) {
  mu = moleculeTarget->mass * gas->mass / (moleculeTarget->mass + gas->mass); 
  cout << "target mass: " << moleculeTarget->mass << " amu" << endl;
  cout << "gas mass: " << gas->mass << " amu" << endl;
  cout << "reduce mass: " << mu << " amu" << endl;
  cout << "charge state: " << moleculeTarget->Q << " e" << endl;
} else if (gas_buffer_flag == 2 || gas_buffer_flag == 3) {
  mu = gas->mass;
  cout << "target mass: " << moleculeTarget->mass << " amu" << endl;
  cout << "gas mass: " << gas->mass << " amu" << endl;
  cout << "reduce mass: " << mu << " amu" << endl;
  cout << "charge state: " << moleculeTarget->Q << " e" << endl;
  Inertia = 0.0;
  for (int i = 0; i < gas->natoms; i++) {
    Inertia += gas->m[i] * pow(gas->z[i],2);
  }
  d_bond = gas->d;
}	

// convert alpha units to kcal/mol
alpha *= ALPHA_TO_KCAL_MOL;

double start_ellipsoid = omp_get_wtime();
if (equipotential_flag == 1) {
  equipotential = new Equipotential(moleculeTarget, gas, polarizability_flag, temperatureTarget, mu, alpha, gas_buffer_flag);
  a = equipotential->a;
  b = equipotential->b;
  c = equipotential->c;
  cout << "ellipsoid axis length: " << a << "  "<< b << "  "<< c << "  Ang" << endl;
} else geometric_ellipsoid();
  
double end_ellipsoid = omp_get_wtime();
cout << "ellipsoid calculation time: " << (end_ellipsoid - start_ellipsoid) << "s" << endl;

bmax = max(max(a,b),c); // maximal impact parameter
 
cout << "maximal impact parameter: " << bmax << " Ang" << endl;

// create the linked cell list
double start_linked_cell = omp_get_wtime();
linkedcell = new LinkedCell(moleculeTarget, a, b, c, lj_cutoff, skin, long_range_flag, long_range_cutoff, coul_cutoff, gas_buffer_flag);
double end_linked_cell = omp_get_wtime(); 
cout << "linked-cell calculation time: " << (end_linked_cell - start_linked_cell) << " s" << endl;

// simulation box length
lx = 0.5*linkedcell->lx;
ly = 0.5*linkedcell->ly;
lz = 0.5*linkedcell->lz;

// loop over iterations
int Niter = nIter;
int Ntraj = nProbe;
double dOmega, Omega, Omega2;
int Nfree, Nscatter, Nlost;

Omega = 0.0;
Omega2 = 0.0;

dOmega = 0.0;
Nscatter = 0;
Nfree = 0;
Nlost = 0;

double start_ccs = omp_get_wtime();

double *dOmega_vec;
double *rnd_vec1, *rnd_vec2, *rnd_vec3, *rnd_vec4, *rnd_vec5;
double *rnd_vec6, *rnd_vec7, *rnd_vec8, *rnd_vec9;
int *Nscatter_vec;
int *Nfree_vec;
int *Nlost_vec;

dOmega_vec = new double [Niter * Ntraj]();
Nscatter_vec = new int [Niter * Ntraj]();
Nfree_vec = new int [Niter * Ntraj]();
Nlost_vec = new int [Niter * Ntraj]();
rnd_vec1 = new double [Niter * Ntraj]();
rnd_vec2 = new double [Niter * Ntraj]();
rnd_vec3 = new double [Niter * Ntraj]();
rnd_vec4 = new double [Niter * Ntraj]();
rnd_vec5 = new double [Niter * Ntraj]();
rnd_vec6 = new double [Niter * Ntraj]();
rnd_vec7 = new double [Niter * Ntraj](); 
rnd_vec8 = new double [Niter * Ntraj]();
rnd_vec9 = new double [Niter * Ntraj](); 


#pragma omp simd
for (int j = 0; j < Niter * Ntraj; j++) {
  rnd_vec1[j] = sqrt(bmax * bmax * mt->getRandomNumber()); // impact parameters vector
  rnd_vec2[j] = mt->getRandomNumber();
  rnd_vec3[j] = mt->getRandomNumber();
  rnd_vec4[j] = mt->getRandomNumber();
  rnd_vec5[j] = mt->getRandomNumber();
  if (gas_buffer_flag == 2 || gas_buffer_flag == 3) {
   rnd_vec6[j] = mt->getRandomNumber();
   rnd_vec7[j] = mt->getRandomNumber();
   rnd_vec8[j] = mt->getRandomNumber();
   rnd_vec9[j] = mt->getRandomNumber();
  }
}

omp_set_num_threads(nthreads);

cout << "*********************************************************" << endl;
cout << "Trajectory calculations " << endl;
cout << "*********************************************************" << endl;
if (gas_buffer_flag == 1 || gas_buffer_flag == 4) { 
  // Hellium: He - atomic 
  #pragma omp parallel for schedule(dynamic)   
  for (int j = 0; j < Niter * Ntraj; j++) {
    bool hit, success;
    double chi;
    Force *force;
    force = new Force(moleculeTarget, linkedcell, lj_cutoff, alpha, coul_cutoff);
    GasBuffer *gasProbe;
    gasProbe = new GasBuffer(gas_buffer_flag);  

    setup(gasProbe, hit, rnd_vec1[j],rnd_vec2[j],rnd_vec3[j],rnd_vec4[j],rnd_vec5[j],rnd_vec6[j],rnd_vec7[j],rnd_vec8[j],rnd_vec9[j]);

    if (hit) {
      run_He(gasProbe, success, chi, dt, force);
      if (success) {
        dOmega_vec[j] = M_PI * (1.0 - cos(chi)) * pow(bmax,2.0);
        Nscatter_vec[j] = 1;
      } else {
        Nlost_vec[j] = 1;
      }
    } else {
      Nfree_vec[j] = 1;
    }
    delete gasProbe;
    delete force;
  }
} else if (gas_buffer_flag == 2) {
  // Nitrogen: N2 - diatomic molecule
  #pragma omp parallel for schedule(dynamic)   
  for (int j = 0; j < Niter * Ntraj; j++) {
    bool hit, success;
    double chi;
    Force *force;
    force = new Force(moleculeTarget, linkedcell, lj_cutoff, alpha, coul_cutoff);
    GasBuffer *gasProbe;
    gasProbe = new GasBuffer(gas_buffer_flag);  

    setup(gasProbe, hit, rnd_vec1[j],rnd_vec2[j],rnd_vec3[j],rnd_vec4[j],rnd_vec5[j],rnd_vec6[j],rnd_vec7[j],rnd_vec8[j],rnd_vec9[j]);

    if (hit) {
      run_N2(gasProbe, success, chi, dt, force);
      if (success) {
        dOmega_vec[j] = M_PI * (1.0 - cos(chi)) * pow(bmax,2.0);
        Nscatter_vec[j] = 1;
      } else {
        Nlost_vec[j] = 1;
      }
    } else {
      Nfree_vec[j] = 1;
    }
    delete gasProbe;
    delete force;
  }
} else if (gas_buffer_flag == 3) {
  // Carbon dioxide: CO2 - linear triatomic molecule
  #pragma omp parallel for schedule(dynamic)   
  for (int j = 0; j < Niter * Ntraj; j++) {
    bool hit, success;
    double chi;
    Force *force;
    force = new Force(moleculeTarget, linkedcell, lj_cutoff, alpha, coul_cutoff);
    GasBuffer *gasProbe;
    gasProbe = new GasBuffer(gas_buffer_flag);  

    setup(gasProbe, hit, rnd_vec1[j],rnd_vec2[j],rnd_vec3[j],rnd_vec4[j],rnd_vec5[j],rnd_vec6[j],rnd_vec7[j],rnd_vec8[j],rnd_vec9[j]);
    
    run_CO2(gasProbe, success, chi, dt, force);
    
    if (hit) {
      run_CO2(gasProbe, success, chi, dt, force);
      if (success) {
        dOmega_vec[j] = M_PI * (1.0 - cos(chi)) * pow(bmax,2.0);
        Nscatter_vec[j] = 1;
      } else {
        Nlost_vec[j] = 1;
      }
    } else {
      Nfree_vec[j] = 1;
    }
    delete gasProbe;
    delete force;
  }
}

int count = 0;
for (int i = 0; i < Niter; i++) {
  dOmega = 0.0;
  Nscatter = 0.0;
  Nfree = 0.0;
  Nlost = 0.0;

  for (int j = 0; j < Ntraj; j++) {
    dOmega += dOmega_vec[count];
    Nscatter += Nscatter_vec[count];
    Nfree += Nfree_vec[count];
    Nlost += Nlost_vec[count];
    count++;
  }
  Omega += 1.0/(float(Nscatter + Nfree))*dOmega;
  Omega2 += pow(1.0/(float(Nscatter + Nfree))*dOmega,2.0);
  printf("Ntraj: %i\n",Ntraj);
  printf("Nfree: %i\n",Nfree);
  printf("Nscatter: %i\n",Nscatter);
  printf("Nlost: %i\n",Nlost);
  printf("omega: %g\n",1.0/(float(Nscatter + Nfree))*dOmega);
}
  
double end_ccs = omp_get_wtime();
cout << "CCS time: " << (end_ccs - start_ccs) << " s" << endl;
  
// average CCS
CCS_ave = Omega/Niter;
double sig2;
sig2 = Omega2/Niter - pow(CCS_ave,2.0);
CCS_err = sqrt(sig2/float(Niter));

cout << "*********************************************************" << endl;
cout << "average value of CCS = " << CCS_ave << " Ang^2" << endl;
cout << "error value of CCS = " << CCS_err << " Ang^2" << endl;

double end = omp_get_wtime();
cout << "Total time: " << (end - start) << " s" << endl;

}

System::~System() {
  delete input;
  delete mt;
  delete moleculeTarget;
  delete gas;
  if(equipotential_flag) delete equipotential;
  delete linkedcell;
}

// Helium gas dynamics
void System::run_He(GasBuffer *gasProbe, bool &success, double &chi, double time_step, Force *force) {
double dt = time_step;
double Ei, Ef, dE, E;
double Ui, Ki, dK, dU;
double Up, Ek;
vector<double> ri(3), rf(3), rcm_i(3);
vector<double> rcm_new(3), rcm_old(3),vcm_new(3),vcm_old(3);
vector<double> vcm_i(3), vcm_f(3);
vector<double> f(3), fi(3);
vector<double> rcm(3), vcm(3);

int maxTries = 8;
int trajTries = 0;
double theta_max = 1.5 * M_PI; // maximal angular displacement
double dtheta = 0.0;
int stepcount;
stepcount = 0;

// copy center mass position and velocities
rcm[0] = gasProbe->x[0];
rcm[1] = gasProbe->y[0];
rcm[2] = gasProbe->z[0];
vcm[0] = gasProbe->vx[0];
vcm[1] = gasProbe->vy[0];
vcm[2] = gasProbe->vz[0];

rcm_i = rcm;
vcm_i = vcm;
rcm_old = rcm_i;
vcm_old = vcm_i;

rcm_new = rcm;

if (force_type == 1) {
  force->lennardjones(gasProbe,0,f,Up);
} else if (force_type == 2) {
  force->lennardjones_LC(gasProbe,0,f,Up);
} else if (force_type == 3) {
  force->lennardjones_induced_dipole(gasProbe,0,f,Up);
} else if (force_type == 4) {
  force->lennardjones_induced_dipole_LC(gasProbe,0,f,Up);
}

fi = f;
Ek = KineticEnergy(mu,vcm);
Ei = Up + Ek;

vector<double> ur(3);
double r;
double tmp, dH;
tmp = 0.0;
int t;
t = 1;

while (trajTries < maxTries && maxTries < 10) {
  
  if (stepcount == 0) {
    // initial conditions;
    vcm = vcm_i;
    rcm = rcm_i;
    rcm_old = rcm_i;
    gasProbe->x[0] = rcm_i[0];
    gasProbe->y[0] = rcm_i[1];
    gasProbe->z[0] = rcm_i[2];
    gasProbe->vx[0] = vcm_i[0];
    gasProbe->vy[0] = vcm_i[1];
    gasProbe->vz[0] = vcm_i[2];
    f = fi;  
    dtheta = 0.0;
    t = 1;
    tmp = 0.0;
  }

  // first-half verlet integration     
  for (int i = 0; i < 3; i++) {
    vcm[i] += 0.5*dt/mu*f[i]*KCALMOLANGAMU_TO_ANGFS2;
    rcm[i] += dt*vcm[i];
  }  
      
  // update gas buffer positions and velocities
  gasProbe->x[0] = rcm[0];
  gasProbe->y[0] = rcm[1];
  gasProbe->z[0] = rcm[2];
  gasProbe->vx[0] = vcm[0];
  gasProbe->vy[0] = vcm[1];
  gasProbe->vz[0] = vcm[2];

  // calculate angular displacement
  rcm_new = rcm;
  dtheta += anglevec(rcm_new,rcm_old);
  rcm_old = rcm_new;

  // calculate force
  if (force_type == 1) {
   force->lennardjones(gasProbe,0,f,Up);
  } else if (force_type == 2) {
   force->lennardjones_LC(gasProbe,0,f,Up);
  } else if (force_type == 3) {
   force->lennardjones_induced_dipole(gasProbe,0,f,Up);
  } else if (force_type == 4) {
   force->lennardjones_induced_dipole_LC(gasProbe,0,f,Up);
  } 

  // second-half verlet integration     
  for (int i = 0; i < 3; i++) {
    vcm[i] += 0.5*dt/mu*f[i]*KCALMOLANGAMU_TO_ANGFS2;
  }

  // update gas buffer velocities
  gasProbe->vx[0] = vcm[0];
  gasProbe->vy[0] = vcm[1];
  gasProbe->vz[0] = vcm[2];
  
  stepcount++;

  if (dtheta > theta_max) {
    trajTries++;
    dt = time_step/(trajTries + 1);
    stepcount = 0;
    continue;
  }
    
  // check if outside the box simulation
  if (abs(rcm[0]) >= lx || abs(rcm[1]) >= ly || abs(rcm[2]) >= lz) {
    // update kinetic energy
    Ek = KineticEnergy(mu,vcm);
    // update total energy
    E = Up + Ek;
    // conservation energy condition
    dE = abs(E-Ei)/Ei*100.0;
    if (dE > 0.5) {  
      trajTries++;
      maxTries++;
      dt = time_step/(trajTries + 1);
      stepcount = 0;
      continue;
    } else {
      vcm_f = vcm;
      chi = anglevec(vcm_i,vcm_f);
      success = true;
      return;
    }
  } else { 
    // check if outside of the ellipsoid
    ur[0] = rcm_new[0]/a;
    ur[1] = rcm_new[1]/b;
    ur[2] = rcm_new[2]/c;
    r = Math::dotProduct(ur,ur);
    if (r > 1.0) {
      // update kinetic energy
      Ek = KineticEnergy(mu,vcm);
      // update total energy
      E = Up + Ek;
      // conservation energy condition
      dE = abs(E-Ei)/Ei*100.0;
      if (dE > 0.5) {  
        trajTries++;
        maxTries++;
        dt = time_step/(trajTries + 1);
        stepcount = 0;
        continue;
      } else {
        vcm_f = vcm;
        chi = anglevec(vcm_i,vcm_f);
        success = true;
        return;
      }
    }
  }
}

success = false;
return;
}

// Nitrogen molecular gas dynamics
void System::run_N2(GasBuffer *gasProbe, bool &success, double &chi, double time_step, Force *force) {
double dt = time_step;
double Ei, Ef, dE, E;
double Ui, Ki, dK, dU;
vector<double> rcm_i(3),rcm_f(3);
vector<double> rcm_new(3),rcm_old(3),vcm_new(3),vcm_old(3);
vector<double> vcm_i(3),vcm_f(3);
vector<double> f(3);
vector<double> rcm(3),vcm(3);

double Up, Ek;
int maxTries = 8;
int trajTries = 0;
double theta_max = 1.5 * M_PI; // maximal angular displacement
double dtheta = 0.0;
int stepcount;
stepcount = 0;

int natoms = gasProbe->natoms;
double x[natoms], y[natoms], z[natoms];
double vx[natoms], vy[natoms], vz[natoms];
double m[natoms];
double M;

// save local gas buffer positions and velocities
for (int iatom = 0; iatom < natoms; iatom++) {
  x[iatom] = gasProbe->x[iatom];
  y[iatom] = gasProbe->y[iatom];
  z[iatom] = gasProbe->z[iatom];
  vx[iatom] = gasProbe->vx[iatom];
  vy[iatom] = gasProbe->vy[iatom];
  vz[iatom] = gasProbe->vz[iatom];
  m[iatom] = gasProbe->m[iatom];
}

M = gasProbe->mass;

// save local intial positions and velocities of buffer gas
double xi[natoms], yi[natoms], zi[natoms];
double vxi[natoms], vyi[natoms], vzi[natoms];
double x_old[natoms], y_old[natoms], z_old[natoms];
double vx_old[natoms], vy_old[natoms], vz_old[natoms];

for (int iatom = 0; iatom < natoms; iatom++) {
  xi[iatom] = x[iatom];
  yi[iatom] = y[iatom];
  zi[iatom] = z[iatom];
  x_old[iatom] = xi[iatom];
  y_old[iatom] = yi[iatom];
  z_old[iatom] = zi[iatom];
  vxi[iatom] = vx[iatom];
  vyi[iatom] = vy[iatom];
  vzi[iatom] = vz[iatom];
  vx_old[iatom] = vxi[iatom];
  vy_old[iatom] = vyi[iatom];
  vz_old[iatom] = vzi[iatom];
}

// copy center mass position
rcm[0] = gasProbe->rcm[0];
rcm[1] = gasProbe->rcm[1];
rcm[2] = gasProbe->rcm[2];

// copy center mass velocity
vcm[0] = gasProbe->vcm[0];
vcm[1] = gasProbe->vcm[1];
vcm[2] = gasProbe->vcm[2];

rcm_i = rcm;
vcm_i = vcm;
rcm_old = rcm_i;
vcm_old = vcm_i;

rcm_new = rcm;

double Up_gas;
double f_gas[natoms][3];
double fi_gas[natoms][3]; // initial force of buffer gas

Up_gas = 0.0;
for (int iatom = 0; iatom < natoms; iatom++) {
  if (force_type == 1) {
    if (iatom == 2) {
      f[0] = 0.0;
      f[1] = 0.0;
      f[2] = 0.0;
      Up = 0.0; 
    } else {
     force->lennardjones(gasProbe,iatom,f,Up);
    } 
  } else if (force_type == 2) {
    if (iatom == 2) {
      f[0] = 0.0;
      f[1] = 0.0;
      f[2] = 0.0;
      Up = 0.0; 
    } else {
      force->lennardjones_LC(gasProbe,iatom,f,Up);
    }
  } else if (force_type == 3) {
    if (iatom == 2) {
      force->coulomb(gasProbe,iatom,f,Up); 
    } else {
      force->lennardjones_coulomb(gasProbe,iatom,f,Up);
    }
  } else if (force_type == 4) {
    if (iatom == 2) {
      force->coulomb_LC(gasProbe,iatom,f,Up); 
    } else {
      force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);
    }
  } else if (force_type == 5) {
    if (iatom == 2) {
      force->coulomb_induced_dipole_iso(gasProbe,iatom,f,Up); 
    } else {
      force->lennardjones_coulomb(gasProbe,iatom,f,Up);
    }
  } else if (force_type == 6) {
    if (iatom == 2) {
      force->coulomb_induced_dipole_iso_LC(gasProbe,iatom,f,Up);
    } else {
      force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);
    }
  }
    
  Up_gas += Up;
  f_gas[iatom][0] = f[0];
  f_gas[iatom][1] = f[1];
  f_gas[iatom][2] = f[2];
  fi_gas[iatom][0] = f[0];
  fi_gas[iatom][1] = f[1];
  fi_gas[iatom][2] = f[2];
}

vector<double> v(3);
Ek = 0.0;
for (int iatom = 0; iatom < natoms; iatom++) {
  v[0] = vx[iatom];
  v[1] = vy[iatom];
  v[2] = vz[iatom];
  Ek += KineticEnergy(m[iatom],v);
}

Ki = Ek;
Ui = Up_gas;

Ei = Up_gas + Ek;

double E_old, E_new;
E_old = Ei;

vector<double> ur(3);
double r;

double tmp, dH;
tmp = 0.0;
int t;
t = 1;

vector<double> ri(3), rj(3), vi(3), vj(3), fi(3), fj(3);
double mi, mj, d2ij;

d2ij = d_bond*d_bond;

double Ek_cm, Ew_rot;
Ek_cm = KineticEnergy(mu,vcm);
Ew_rot = Ek - Ek_cm;

while (trajTries < maxTries && maxTries < 10) {
  
  if (stepcount == 0) {
    // initial conditions;
    vcm = vcm_i;
    rcm = rcm_i;
    rcm_old = rcm_i;
    gasProbe->rcm[0] = rcm[0];
    gasProbe->rcm[1] = rcm[1];
    gasProbe->rcm[2] = rcm[2];
    gasProbe->vcm[0] = vcm[0];
    gasProbe->vcm[1] = vcm[1];
    gasProbe->vcm[2] = vcm[2];
    for (int iatom = 0; iatom < natoms; iatom++) {
      gasProbe->x[iatom] = xi[iatom];
      gasProbe->y[iatom] = yi[iatom];
      gasProbe->z[iatom] = zi[iatom];
      gasProbe->vx[iatom] = vxi[iatom];
      gasProbe->vy[iatom] = vyi[iatom];
      gasProbe->vz[iatom] = vzi[iatom];
      x[iatom] = xi[iatom];
      y[iatom] = yi[iatom];
      z[iatom] = zi[iatom];
      vx[iatom] = vxi[iatom];
      vy[iatom] = vyi[iatom];
      vz[iatom] = vzi[iatom];
      for (int i = 0; i < 3; i++) {
        f_gas[iatom][i] = fi_gas[iatom][i];
      }
    }
    dtheta = 0.0;
    t = 1;
    tmp = 0.0;
  }

  // first-half verlet integration     
  // N2 diatomic molecule N-N
  // nitrogen 1
  ri[0] = x[0];
  ri[1] = y[0];
  ri[2] = z[0];
  vi[0] = vx[0];
  vi[1] = vy[0];
  vi[2] = vz[0];
  mi = m[0];
  // nitrogen 2
  rj[0] = x[1];
  rj[1] = y[1];
  rj[2] = z[1];
  vj[0] = vx[1];
  vj[1] = vy[1];
  vj[2] = vz[1];
  mj = m[1];

  for (int i = 0; i < 3; i++) {
    fi[i] = f_gas[0][i] + 0.5*f_gas[2][i]; // nitrogen 1 + dummy atom
    fj[i] = f_gas[1][i] + 0.5*f_gas[2][i]; // nitrogen 2 + dummy atom
  }   
    
  first_half_verlet_constrained(ri, rj, vi, vj, fi, fj, mi, mj, dt, d2ij);

  // update positions and velocties
  // nitrogen 1
  x[0] = ri[0];
  y[0] = ri[1];
  z[0] = ri[2];
  vx[0] = vi[0];
  vy[0] = vi[1];
  vz[0] = vi[2];
  // nitrogen 2
  x[1] = rj[0];
  y[1] = rj[1];
  z[1] = rj[2];
  vx[1] = vj[0];
  vy[1] = vj[1];
  vz[1] = vj[2];
  // update center mass position and velocity
  for (int i = 0; i < 3; i++) {
    rcm[i] = 0.5*(ri[i] + rj[i]);
    vcm[i] = 0.5*(vi[i] + vj[i]);      
  }
  // dummy
  x[2] = rcm[0];
  y[2] = rcm[1];
  z[2] = rcm[2];
  
  // update gas buffer positions and velocities
  for (int iatom = 0; iatom < natoms; iatom++) {
    gasProbe->x[iatom] = x[iatom];
    gasProbe->y[iatom] = y[iatom];
    gasProbe->z[iatom] = z[iatom];
    gasProbe->vx[iatom] = vx[iatom];
    gasProbe->vy[iatom] = vy[iatom];
    gasProbe->vz[iatom] = vz[iatom]; 
  }

  // calculate angular displacement
  rcm_new = rcm;
  dtheta += anglevec(rcm_new,rcm_old);
  rcm_old = rcm_new;

  // calculate force
  Up_gas = 0.0;
  for (int iatom = 0; iatom < natoms; iatom++) {
    if (force_type == 1) {
      if (iatom == 2) {
        f[0] = 0.0;
        f[1] = 0.0;
        f[2] = 0.0;
        Up = 0.0; 
      } else {
        force->lennardjones(gasProbe,iatom,f,Up);
      } 
    } else if (force_type == 2) {
      if (iatom == 2) {
        f[0] = 0.0;
        f[1] = 0.0;
        f[2] = 0.0;
        Up = 0.0; 
      } else {
        force->lennardjones_LC(gasProbe,iatom,f,Up);
      }
    } else if (force_type == 3) {
      if (iatom == 2) {
        force->coulomb(gasProbe,iatom,f,Up); 
      } else {
        force->lennardjones_coulomb(gasProbe,iatom,f,Up);
      }
    } else if (force_type == 4) {
      if (iatom == 2) {
        force->coulomb_LC(gasProbe,iatom,f,Up); 
      } else {
        force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);
      }
    } else if (force_type == 5) {
      if (iatom == 2) {
        force->coulomb_induced_dipole_iso(gasProbe,iatom,f,Up);
      } else {
        force->lennardjones_coulomb(gasProbe,iatom,f,Up);
      }
    } else if (force_type == 6) {
      if (iatom == 2) {
        force->coulomb_induced_dipole_iso_LC(gasProbe,iatom,f,Up);
      } else {
        force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);
      }
    }  
    Up_gas += Up;
    f_gas[iatom][0] = f[0];
    f_gas[iatom][1] = f[1];
    f_gas[iatom][2] = f[2];
  } 
  
  // second-half verlet integration     
  // nitrogen 1
  ri[0] = x[0];
  ri[1] = y[0];
  ri[2] = z[0];
  vi[0] = vx[0];
  vi[1] = vy[0];
  vi[2] = vz[0];
  mi = m[0];
  // nitrogen 2
  rj[0] = x[1];
  rj[1] = y[1];
  rj[2] = z[1];
  vj[0] = vx[1];
  vj[1] = vy[1];
  vj[2] = vz[1];
  mj = m[1];

  for (int i = 0; i < 3; i++) {
    fi[i] = f_gas[0][i] + 0.5*f_gas[2][i]; // nitrogen 1 + dummy atom
    fj[i] = f_gas[1][i] + 0.5*f_gas[2][i]; // nitrogen 2 + dummy atom
  }   
    
  second_half_verlet_constrained(ri, rj, vi, vj, fi, fj, mi, mj, dt);
  // update velocities
  vx[0] = vi[0];
  vy[0] = vi[1];
  vz[0] = vi[2];
  vx[1] = vj[0];
  vy[1] = vj[1];
  vz[1] = vj[2];

  // update center mass position and velocity
  for (int i = 0; i < 3; i++) {
    vcm[i] = 0.5*(vi[i] + vj[i]);      
  }

  // update gas buffer velocities
  for (int iatom = 0; iatom < natoms-1; iatom++) {
    gasProbe->vx[iatom] = vx[iatom];
    gasProbe->vy[iatom] = vy[iatom];
    gasProbe->vz[iatom] = vz[iatom]; 
  }
  
  stepcount++;
   
  if (dtheta > theta_max) {
    trajTries++;
    dt = time_step/(trajTries + 1);
    stepcount = 0;
    continue;
  }
  
  // check if outside the box simulation
  if (abs(rcm[0]) >= lx || abs(rcm[1]) >= ly || abs(rcm[2]) >= lz) {
    Ek = 0.0;
    for (int iatom = 0; iatom < natoms-1; iatom++) {
      v[0] = vx[iatom];
      v[1] = vy[iatom];
      v[2] = vz[iatom];
      Ek += KineticEnergy(m[iatom],v);
    }
    E = Up_gas + Ek;
    dE = abs((E-Ei)/Ei)*100.0;
    if (dE > 1.0) {
      trajTries++;
      maxTries++;
      dt = time_step/(trajTries + 1);
      stepcount = 0;
      continue;
    } else { 
      vcm_f = vcm;
      chi = anglevec(vcm_i,vcm_f);
      success = true;
      return;
    }
  } else { 
    // check if outside of the ellipsoid
    ur[0] = rcm_new[0]/a;
    ur[1] = rcm_new[1]/b;
    ur[2] = rcm_new[2]/c;
    r = Math::dotProduct(ur,ur);
    if (r > 1.0) {
      Ek = 0.0;
      for (int iatom = 0; iatom < natoms-1; iatom++) {
        v[0] = vx[iatom];
        v[1] = vy[iatom];
        v[2] = vz[iatom];
        Ek += KineticEnergy(m[iatom],v);
      }
      E = Up_gas + Ek;
      dE = abs((E-Ei)/Ei)*100.0;
      if (dE > 1.0) {
        trajTries++;
        maxTries++;
        dt = time_step/(trajTries + 1);
        stepcount = 0;
        continue;
      } else {
        vcm_f = vcm;
        chi = anglevec(vcm_i,vcm_f);
        success = true;
        return;
      }
    }
  }
}

success = false;
return;
}

// Carbon dioxide molecular gas dynamics
void System::run_CO2(GasBuffer *gasProbe, bool &success, double &chi, double time_step, Force *force) {
double dt = time_step;
double Ei, Ef, dE, E, Wc, W;
double Ui, Ki, dK, dU;
vector<double> rcm_i(3),rcm_f(3);
vector<double> rcm_new(3),rcm_old(3),vcm_new(3),vcm_old(3);
vector<double> vcm_i(3),vcm_f(3);
vector<double> f(3);
vector<double> rcm(3),vcm(3);

double Up, Ek;
int maxTries = 8;
int trajTries = 0;
double theta_max = 1.5 * M_PI; // maximal angular displacement
double dtheta = 0.0;
int stepcount;
stepcount = 0;

int natoms = gasProbe->natoms;
double x[natoms], y[natoms], z[natoms];
double vx[natoms], vy[natoms], vz[natoms];
double m[natoms];
double M;

// save local gas buffer positions and velocities
for (int iatom = 0; iatom < natoms; iatom++) {
  x[iatom] = gasProbe->x[iatom];
  y[iatom] = gasProbe->y[iatom];
  z[iatom] = gasProbe->z[iatom];
  vx[iatom] = gasProbe->vx[iatom];
  vy[iatom] = gasProbe->vy[iatom];
  vz[iatom] = gasProbe->vz[iatom];
  m[iatom] = gasProbe->m[iatom];
}

M = gasProbe->mass;

// save local intial positions and velocities of buffer gas
double xi[natoms], yi[natoms], zi[natoms];
double vxi[natoms], vyi[natoms], vzi[natoms];
double x_old[natoms], y_old[natoms], z_old[natoms];
double vx_old[natoms], vy_old[natoms], vz_old[natoms];

for (int iatom = 0; iatom < natoms; iatom++) {
  xi[iatom] = x[iatom];
  yi[iatom] = y[iatom];
  zi[iatom] = z[iatom];
  x_old[iatom] = xi[iatom];
  y_old[iatom] = yi[iatom];
  z_old[iatom] = zi[iatom];
  vxi[iatom] = vx[iatom];
  vyi[iatom] = vy[iatom];
  vzi[iatom] = vz[iatom];
  vx_old[iatom] = vxi[iatom];
  vy_old[iatom] = vyi[iatom];
  vz_old[iatom] = vzi[iatom];
}

// copy center mass position
rcm[0] = gasProbe->rcm[0];
rcm[1] = gasProbe->rcm[1];
rcm[2] = gasProbe->rcm[2];

// copy center mass velocity
vcm[0] = gasProbe->vcm[1];
vcm[1] = gasProbe->vcm[1];
vcm[2] = gasProbe->vcm[2];

rcm_i = rcm;
vcm_i = vcm;
rcm_old = rcm_i;
vcm_old = vcm_i;

rcm_new = rcm;

double Up_gas;
double f_gas[natoms][3];
double fi_gas[natoms][3]; // initial force of buffer gas

Up_gas = 0.0;
for (int iatom = 0; iatom < natoms; iatom++) {
  if (force_type == 1) {
    if (iatom == 1) { 
      force->lennardjones_CO2(gasProbe,iatom,f,Up);
    } else {
      force->lennardjones(gasProbe,iatom,f,Up);
    }
  } else if (force_type == 2) {
    if (iatom == 1) { 
      force->lennardjones_LC_CO2(gasProbe,iatom,f,Up);
    } else {
      force->lennardjones_LC(gasProbe,iatom,f,Up);
    }    
  } else if (force_type == 3) {
    if (iatom == 1) { 
      force->lennardjones_coulomb_CO2(gasProbe,iatom,f,Up);    
    } else {
      force->lennardjones_coulomb(gasProbe,iatom,f,Up);    
    }  
  } else if (force_type == 4) {
    if (iatom == 1) { 
      force->lennardjones_coulomb_LC_CO2(gasProbe,iatom,f,Up);    
    } else {
      force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);    
    }
  } else if (force_type == 5) {
    if (iatom == 1) {
      force->lennardjones_coulomb_induced_dipole_iso_CO2(gasProbe,iatom,f,Up);
    } else {
      force->lennardjones_coulomb(gasProbe,iatom,f,Up);
    }
  } else if (force_type == 6) {
    if (iatom == 1) {
      force->lennardjones_coulomb_induced_dipole_iso_LC_CO2(gasProbe,iatom,f,Up);
    } else {
      force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);
    }
  }
  
  Up_gas += Up;
  f_gas[iatom][0] = f[0];
  f_gas[iatom][1] = f[1];
  f_gas[iatom][2] = f[2];
  fi_gas[iatom][0] = f[0];
  fi_gas[iatom][1] = f[1];
  fi_gas[iatom][2] = f[2];
}

vector<double> v(3);
Ek = 0.0;
for (int iatom = 0; iatom < natoms; iatom++) {
  v[0] = vx[iatom];
  v[1] = vy[iatom];
  v[2] = vz[iatom];
  Ek += KineticEnergy(m[iatom],v);
}

Ki = Ek;
Ui = Up_gas;

Ei = Up_gas + Ek;

double E_old, E_new;
E_old = Ei;

vector<double> ur(3);
double r;

double tmp, dH;
tmp = 0.0;
int t;
t = 1;

vector<double> ri(3), rj(3), vi(3), vj(3), fi(3), fj(3);
double mi, mj, m2M, d2ij;

d2ij = d_bond*d_bond;

double Ek_cm, Ew_rot;
Ek_cm = KineticEnergy(mu,vcm);
Ew_rot = Ek - Ek_cm;

while (trajTries < maxTries && maxTries < 10) {
  
  if (stepcount == 0) {
    // initial conditions;
    vcm = vcm_i;
    rcm = rcm_i;
    rcm_old = rcm_i;
    gasProbe->rcm[0] = rcm[0];
    gasProbe->rcm[1] = rcm[1];
    gasProbe->rcm[2] = rcm[2];
    gasProbe->vcm[0] = vcm[0];
    gasProbe->vcm[1] = vcm[1];
    gasProbe->vcm[2] = vcm[2];
    for (int iatom = 0; iatom < natoms; iatom++) {
      gasProbe->x[iatom] = xi[iatom];
      gasProbe->y[iatom] = yi[iatom];
      gasProbe->z[iatom] = zi[iatom];
      gasProbe->vx[iatom] = vxi[iatom];
      gasProbe->vy[iatom] = vyi[iatom];
      gasProbe->vz[iatom] = vzi[iatom];
      x[iatom] = xi[iatom];
      y[iatom] = yi[iatom];
      z[iatom] = zi[iatom];
      vx[iatom] = vxi[iatom];
      vy[iatom] = vyi[iatom];
      vz[iatom] = vzi[iatom];
      for (int i = 0; i < 3; i++) {
        f_gas[iatom][i] = fi_gas[iatom][i];
      }
    }
    dtheta = 0.0;
    t = 1;
    tmp = 0.0;
  }

  // first-half verlet integration     
  // CO2 trilinear molecule O-C-O
  // G. Ciccotti et al. Molecular dynamics of rigid systems 
  // in cartesian coordinates A general formulation
  // MOLECULAR PHYSICS, 1982, VOL. 47, No. 6, 1253-1264
  m2M = m[1]/M;
  for (int i = 0; i < 3; i++) {
    fi[i] = (1.0 - 0.5*m2M)*f_gas[0][i] + m[0]/M*f_gas[1][i] - 0.5*m2M*f_gas[2][i];
    fj[i] = (1.0 - 0.5*m2M)*f_gas[2][i] + m[0]/M*f_gas[1][i] - 0.5*m2M*f_gas[0][i];
  }  
  // oxygen 1
  ri[0] = x[0];
  ri[1] = y[0];
  ri[2] = z[0];
  vi[0] = vx[0];
  vi[1] = vy[0];
  vi[2] = vz[0];
  mi = m[0];
  // oxyegn 2
  rj[0] = x[2];
  rj[1] = y[2];
  rj[2] = z[2];
  vj[0] = vx[2];
  vj[1] = vy[2];
  vj[2] = vz[2];
  mj = m[2];

  first_half_verlet_constrained(ri, rj, vi, vj, fi, fj, mi, mj, dt, d2ij);
  // update positions and velocties
  // oxygen 1
  x[0] = ri[0];
  y[0] = ri[1];
  z[0] = ri[2];
  vx[0] = vi[0];
  vy[0] = vi[1];
  vz[0] = vi[2];
  // oxygen 2
  x[2] = rj[0];
  y[2] = rj[1];
  z[2] = rj[2];
  vx[2] = vj[0];
  vy[2] = vj[1];
  vz[2] = vj[2];
  // update center mass position and velocity
  for (int i = 0; i < 3; i++) {
    rcm[i] = 0.5*(ri[i]+rj[i]);
    vcm[i] = 0.5*(vi[i]+vj[i]);
  }
  // carbon
  x[1] = rcm[0];
  y[1] = rcm[1];
  z[1] = rcm[2];
  vx[1] = vcm[0];
  vy[1] = vcm[1];
  vz[1] = vcm[2];

  // update gas buffer positions and velocities
  for (int iatom = 0; iatom < natoms; iatom++) {
    gasProbe->x[iatom] = x[iatom];
    gasProbe->y[iatom] = y[iatom];
    gasProbe->z[iatom] = z[iatom];
    gasProbe->vx[iatom] = vx[iatom];
    gasProbe->vy[iatom] = vy[iatom];
    gasProbe->vz[iatom] = vz[iatom]; 
  }

  // calculate angular displacement
  rcm_new = rcm;
  dtheta += anglevec(rcm_new,rcm_old);
  rcm_old = rcm_new;
  
  // calculate force
  Up_gas = 0.0;
  
  for (int iatom = 0; iatom < natoms; iatom++) {
    if (force_type == 1) {
      if (iatom == 1) { 
        force->lennardjones_CO2(gasProbe,iatom,f,Up);
      } else {
        force->lennardjones(gasProbe,iatom,f,Up);
      }
    } else if (force_type == 2) {
      if (iatom == 1) { 
        force->lennardjones_LC_CO2(gasProbe,iatom,f,Up);
      } else {
        force->lennardjones_LC(gasProbe,iatom,f,Up);
      }    
    } else if (force_type == 3) {
      if (iatom == 1) { 
        force->lennardjones_coulomb_CO2(gasProbe,iatom,f,Up);    
      } else {
        force->lennardjones_coulomb(gasProbe,iatom,f,Up);    
      }  
    } else if (force_type == 4) {
      if (iatom == 1) { 
        force->lennardjones_coulomb_LC_CO2(gasProbe,iatom,f,Up);    
      } else {
        force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);    
      }
    } else if (force_type == 5) {
      if (iatom == 1) {
        force->lennardjones_coulomb_induced_dipole_iso_CO2(gasProbe,iatom,f,Up);
      } else {
        force->lennardjones_coulomb(gasProbe,iatom,f,Up);
      }
    } else if (force_type == 6) {
      if (iatom == 1) {
        force->lennardjones_coulomb_induced_dipole_iso_LC_CO2(gasProbe,iatom,f,Up);
      } else {
        force->lennardjones_coulomb_LC(gasProbe,iatom,f,Up);
      }
    }
  
    Up_gas += Up;
    f_gas[iatom][0] = f[0];
    f_gas[iatom][1] = f[1];
    f_gas[iatom][2] = f[2];
  } 

  // second-half verlet integration     
  m2M = m[1]/M;
  for (int i = 0; i < 3; i++) {
    fi[i] = (1.0 - 0.5*m2M)*f_gas[0][i] + m[0]/M*f_gas[1][i] - 0.5*m2M*f_gas[2][i];
    fj[i] = (1.0 - 0.5*m2M)*f_gas[2][i] + m[0]/M*f_gas[1][i] - 0.5*m2M*f_gas[0][i];
  }
  // oxygen 1
  ri[0] = x[0];
  ri[1] = y[0];
  ri[2] = z[0];
  vi[0] = vx[0];
  vi[1] = vy[0];
  vi[2] = vz[0];
  mi = m[0];
  // oxyegn 2
  rj[0] = x[2];
  rj[1] = y[2];
  rj[2] = z[2];
  vj[0] = vx[2];
  vj[1] = vy[2];
  vj[2] = vz[2];
  mj = m[2];

  second_half_verlet_constrained(ri, rj, vi, vj, fi, fj, mi, mj, dt);
  // oxygen 1
  vx[0] = vi[0];
  vy[0] = vi[1];
  vz[0] = vi[2];
  // oxygen 2
  vx[2] = vj[0];
  vy[2] = vj[1];
  vz[2] = vj[2];
  // update center mass position and velocity
  for (int i = 0; i < 3; i++) {
    vcm[i] = 0.5*(vi[i]+vj[i]);
  }
  // carbon
  vx[1] = vcm[0];
  vy[1] = vcm[1];
  vz[1] = vcm[2];

  // update gas buffer velocities
  for (int iatom = 0; iatom < natoms; iatom++) {
    gasProbe->vx[iatom] = vx[iatom];
    gasProbe->vy[iatom] = vy[iatom];
    gasProbe->vz[iatom] = vz[iatom]; 
  }
  
  stepcount++;

  if (dtheta > theta_max) {
    trajTries++;
    dt = time_step/(trajTries + 1);
    stepcount = 0;
    continue;
  } 

  // check if outside the box simulation
  if (abs(rcm[0]) >= lx || abs(rcm[1]) >= ly || abs(rcm[2]) >= lz) {
    Ek = 0.0;
    for (int iatom = 0; iatom < natoms; iatom++) {
      v[0] = vx[iatom];
      v[1] = vy[iatom];
      v[2] = vz[iatom];
      Ek += KineticEnergy(m[iatom],v);
    }
    E = Up_gas + Ek;
    dE = abs((E-Ei)/Ei)*100.0;
    if (dE > 5.0) {
      trajTries++;
      maxTries++;
      dt = time_step/(trajTries + 1);
      stepcount = 0;
      continue;
    } else {
      vcm_f = vcm;
      chi = anglevec(vcm_i,vcm_f);
      success = true;
      return;
    }
  } else { 
    // check if outside of the ellipsoid
    ur[0] = rcm_new[0]/a;
    ur[1] = rcm_new[1]/b;
    ur[2] = rcm_new[2]/c;
    r = Math::dotProduct(ur,ur);
    if (r > 1.0) {
      Ek = 0.0;
      for (int iatom = 0; iatom < natoms; iatom++) {
        v[0] = vx[iatom];
        v[1] = vy[iatom];
        v[2] = vz[iatom];
        Ek += KineticEnergy(m[iatom],v);
      }
      E = Up_gas + Ek;
      dE = abs((E-Ei)/Ei)*100.0;      
      if (dE > 5.0) {
        trajTries++;
        maxTries++;
        dt = time_step/(trajTries + 1);
        stepcount = 0;
        continue;
      } else {
        vcm_f = vcm;
        chi = anglevec(vcm_i,vcm_f);
        success = true;
        return;
      }
    }
  }
}

success = false;
return;
}

/* 
 * Inicial position and velocity of buffer gas on ellipsoid surface
 */
void System::setup(GasBuffer *gasProbe, bool &hit, double rndVal1, double rndVal2, double rndVal3, double rndVal4, 
double rndVal5, double rndVal6, double rndVal7, double rndVal8, double rndVal9) {
double xProbe, yProbe, zProbe, gamma, phi, theta, bi;
vector<double> rcm(3);
vector<double> vcm(3);
vector<double> angles(3);
// impact parameter
bi = rndVal1;

// define intial position of the molecule probe
rcm[0] = bi;
rcm[1] = 0.0;
rcm[2] = bmax;

// define intial velocity of the molecule probe
double vi;
vi = velGenerator(mu, temperatureTarget, rndVal2); // velocity distribution
vcm[0] = 0.0;
vcm[1] = 0.0;
vcm[2] = -vi; // velocity in the z direction

// rotational angles
gamma = 2.0 * M_PI * rndVal3; // angle of 0 to 2 Pi
phi = asin(2.0 * rndVal4 - 1.0) + 0.5*M_PI; // angle of -Pi/2 to Pi/2
theta = 2.0 * M_PI * rndVal5; // angle of 0 to 2 Pi
angles[0] = gamma;
angles[1] = phi;
angles[2] = theta;

// rotate position and velocity of molecule probe
rotate(rcm, vcm, angles);

// hit the ellipsoid
vector<double> ur(3), uv(3), cross(3);
double delta, lambda, l1, l2;

ur[0] = rcm[0]/a;
ur[1] = rcm[1]/b;
ur[2] = rcm[2]/c;

uv[0] = vcm[0]/a;
uv[1] = vcm[1]/b;
uv[2] = vcm[2]/c;

cross = Math::crossProduct(ur,uv);
delta = -Math::dotProduct(cross,cross) + Math::dotProduct(uv,uv);

if (delta < 0) {
  hit = false;
  return;
}

lambda = - Math::dotProduct(ur,uv)/Math::dotProduct(uv,uv) - sqrt(delta)/Math::dotProduct(uv,uv);

vector<double> d(3);

d[0] = lambda*vcm[0];
d[1] = lambda*vcm[1];
d[2] = lambda*vcm[2];

rcm[0] = rcm[0] + d[0];
rcm[1] = rcm[1] + d[1];
rcm[2] = rcm[2] + d[2];

// position and velocity of center mass gas buffer
gasProbe->rcm[0] = rcm[0];
gasProbe->rcm[1] = rcm[1];
gasProbe->rcm[2] = rcm[2];
gasProbe->vcm[0] = vcm[0];
gasProbe->vcm[1] = vcm[1];
gasProbe->vcm[2] = vcm[2];

if (gas_buffer_flag == 1 || gas_buffer_flag == 4) {
  gasProbe->x[0] = rcm[0];
  gasProbe->y[0] = rcm[1];
  gasProbe->z[0] = rcm[2];
  gasProbe->vx[0] = vcm[0];
  gasProbe->vy[0] = vcm[1];
  gasProbe->vz[0] = vcm[2];
} else if (gas_buffer_flag == 2 || gas_buffer_flag == 3) {
  // angular position	
  theta = asin(2.0 * rndVal6 - 1.0) + 0.5*M_PI; // angle of -Pi/2 to Pi/2
  phi = 2.0*M_PI * rndVal7; // angle of 0 to 2Pi
  // angular velocity
  double psi, omega;
  psi = 2.0*M_PI * rndVal8; // angle of 0 to 2Pi
  omega = sqrt(2.0*BOLTZMANN_K * temperatureTarget / Inertia * log(1./(1.-rndVal9))) * OMEGA_TO_FS_INV; // angular velocity distribution
  vector<double> rgas(3), wgas(3), vgas(3);
  wgas[0] = omega * cos(psi);
  wgas[1] = omega * sin(psi);
  wgas[2] = 0.0;  
   
  rotate_gas(wgas, theta, phi);
  int natoms = gas->natoms;
  for (int i = 0; i < natoms; i++) {
    rgas[0] = gasProbe->x[i];
    rgas[1] = gasProbe->y[i];
    rgas[2] = gasProbe->z[i];
    rotate_gas(rgas, theta, phi);
    gasProbe->x[i] = rgas[0];
    gasProbe->y[i] = rgas[1];
    gasProbe->z[i] = rgas[2];
    vgas = Math::crossProduct(wgas,rgas);
    gasProbe->vx[i] = vgas[0];
    gasProbe->vy[i] = vgas[1];
    gasProbe->vz[i] = vgas[2];  
  }

  // translate to ellipsoid surface
  for (int i = 0; i < natoms; i++) {
    gasProbe->x[i] += rcm[0];
    gasProbe->y[i] += rcm[1];
    gasProbe->z[i] += rcm[2];
    gasProbe->vx[i] += vcm[0];
    gasProbe->vy[i] += vcm[1];
    gasProbe->vz[i] += vcm[2];
  }     
} 

hit = true;
return;
}

/* 
 * velocity distribution of buffer gas
 */
double System::velDistr(double v, double m, double temperature) {
// convert mass from au to kg
m = m * AMU_TO_KG;

// convert v from Ang/fs to m/s
v = v * ANG_TO_M / FS_TO_S;

return pow(v, 5.0) * pow((m / (2.0 * BOLTZMANN_K * temperature)), 3.0) *
  exp(-(m * pow(v, 2.0)) / (2.0 * BOLTZMANN_K * temperature));
}

/* 
 * velocity distribution of buffer gas
 */
double System::velGenerator(double m, double temperature, double sd) {
//
double v, dv;
double sumProbability;

// initialize v, dv and m
v = 0.0; // units A/fs
dv = 1.0E-5;

// initial probability sum
sumProbability = velDistr(v, m, temperature);

while (sumProbability < sd) {
  v += dv;
  sumProbability += velDistr(v, m, temperature);
}

return v;
}

/* 
 * kinetic energy of buffer gas
 */
double System::KineticEnergy(double m, vector<double> v) {
double Ek;
//convert from amu*ang^2*fs^-2 to kcal/mol
return 0.5*m*Math::dotProduct(v,v)*amuAngfs2_to_KCAL_MOL;
}

double System::anglevec(vector<double> vi, vector<double> vf) {
double t1, t2, t3, argument;

t1 = Math::dotProduct(vi, vf);
t2 = Math::vecModulus(vi);
t3 = Math::vecModulus(vf);

argument = t1 / (t2 * t3);

if (fabs(argument) > 1.0)
return 0.0;

return acos(argument);
}
/* 
 * rotation center masss over ellipsoid surface
 */
void System::rotate(vector<double> &r, vector<double> &v, vector<double> angles) {
vector<double> ri(3,0);
vector<double> vi(3,0);

ri = r;
vi = v;
// Rz	
r[0] = ri[0] * cos(angles[0]) - ri[1] * sin(angles[0]);
r[1] = ri[0] * sin(angles[0]) + ri[1] * cos(angles[0]);
r[2] = ri[2];

v[0] = vi[0] * cos(angles[0]) - vi[1] * sin(angles[0]);
v[1] = vi[0] * sin(angles[0]) + vi[1] * cos(angles[0]);
v[2] = vi[2];

ri = r;
vi = v;
// Ry
r[0] = ri[0] * cos(angles[1]) + ri[2] * sin(angles[1]);
r[1] = ri[1];
r[2] = ri[0] * -sin(angles[1]) + ri[2] * cos(angles[1]);

v[0] = vi[0] * cos(angles[1]) + vi[2] * sin(angles[1]);
v[1] = vi[1];
v[2] = vi[0] * -sin(angles[1]) + vi[2] * cos(angles[1]);

ri = r;
vi = v;
// Rz	
r[0] = ri[0] * cos(angles[2]) - ri[1] * sin(angles[2]);
r[1] = ri[0] * sin(angles[2]) + ri[1] * cos(angles[2]);
r[2] = ri[2];

v[0] = vi[0] * cos(angles[2]) - vi[1] * sin(angles[2]);
v[1] = vi[0] * sin(angles[2]) + vi[1] * cos(angles[2]);
v[2] = vi[2];
}

/* 
 * rotation gas buffer
 */
void System::rotate_gas(vector<double> &v, double theta, double phi) {
vector<double> vi(3);

vi = v;
// Ry
v[0] = vi[0] * cos(theta) + vi[2] * sin(theta);
v[1] = vi[1];
v[2] = vi[0] * -sin(theta) + vi[2] * cos(theta);

vi = v;
// Rz	
v[0] = vi[0] * cos(phi) - vi[1] * sin(phi);
v[1] = vi[0] * sin(phi) + vi[1] * cos(phi);
v[2] = vi[2];
}

void System::geometric_ellipsoid() {
double x, y, z, d;
double xi, yi, zi, ri, ux, uy, uz;
double maxX, maxY, maxZ;
double ellipsoid;
double radius;
maxX = abs(moleculeTarget->maxX);
maxY = abs(moleculeTarget->maxY);
maxZ = abs(moleculeTarget->maxZ);
d = gas->d; 

if (long_range_flag == 1) {
  radius = 2.0*coul_cutoff;
} else {
  radius = 2.0*lj_cutoff;
}
cout << "*********************************************************" << endl;
cout << "Geometric Ellipsoid: " << endl;
cout << "*********************************************************" << endl;
cout << "maximal distances: " << maxX << "  "<< maxY << "  "<< maxZ << "  Ang" << endl;

a = maxX + radius + skin + d;
b = maxY + radius + skin + d;
c = maxZ + radius + skin + d;

cout << "Initial axis length: " << a << "  "<< b << "  "<< c << "  Ang" << endl;

double delta = 0.01;

bool outside;

for (int i = 0; i < moleculeTarget->natoms; i++) {
  xi = moleculeTarget->x[i];
  yi = moleculeTarget->y[i];
  zi = moleculeTarget->z[i];
  ri = sqrt(xi*xi + yi*yi + zi*zi);
  ux = xi/ri;
  uy = yi/ri;
  uz = zi/ri;
  x = xi + (radius+d)*ux;
  y = yi + (radius+d)*uy;
  z = zi + (radius+d)*uz;

  do {
    ellipsoid = pow(x/a,2) + pow(y/b,2) + pow(z/c,2);

    if (ellipsoid > 1.0) {
      a += delta;
      b += delta;
      c += delta;
      outside = true;
    } else {
      outside = false;
    }
  } while(outside);
}

cout << "Ellipsoid axis length: " << a << "  "<< b << "  "<< c << "  Ang" << endl;
}
	
void System::first_half_verlet_constrained(vector<double> &ri, vector<double> &rj, vector<double> &vi, vector<double> &vj,
vector<double> fi, vector<double> fj, double mi, double mj, double dt, double d2ij) {
double muij, delta, g;
int it;
vector<double> rci(3), rcj(3), vci(3), vcj(3), rcij(3), rij(3), gij(3); 

muij = mi*mj/(mi+mj);

for (int i = 0; i < 3; i++) {
  vci[i] = vi[i] + 0.5*dt/mi*fi[i]*KCALMOLANGAMU_TO_ANGFS2; // v'1(dt/2)
  vcj[i] = vj[i] + 0.5*dt/mj*fj[i]*KCALMOLANGAMU_TO_ANGFS2; // v'2(dt/2)
  rci[i] = ri[i] + dt*vci[i];                               // r'1(dt) 
  rcj[i] = rj[i] + dt*vcj[i];                               // r'2(dt)
  rij[i] = ri[i] - rj[i];                                   // r12(0)  
  rcij[i] = rci[i] - rcj[i];                                // r'12(dt)  
}  

double rc2ij,r2ij,rcijrij;
r2ij = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
rc2ij = rcij[0]*rcij[0] + rcij[1]*rcij[1] + rcij[2]*rcij[2];
rcijrij = rcij[0]*rij[0] + rcij[1]*rij[1] + rcij[2]*rij[2];

double determ;
determ = rcijrij*rcijrij - r2ij*(rc2ij - d2ij);
if (determ < 0.0) {
  determ = 0.0;
} 

double lambda;
lambda = (rcijrij - sqrt(determ))/r2ij; // adimensional parameter
// update position and velocities after apply first constrained

for (int i = 0; i < 3; i++) {
  vi[i] = vci[i] - muij/(mi*dt)*lambda*rij[i]; // v1(dt/2)
  vj[i] = vcj[i] + muij/(mj*dt)*lambda*rij[i]; // v2(dt/2)
  ri[i] = rci[i] - muij/mi*lambda*rij[i];      // r1(dt)
  rj[i] = rcj[i] + muij/mi*lambda*rij[i];      // r2(dt)
}

}

void System::second_half_verlet_constrained(vector<double> ri, vector<double> rj, vector<double> &vi, vector<double> &vj,
 vector<double> fi, vector<double> fj, double mi, double mj, double dt) {
double muij, delta, g;
int it;
vector<double> rij(3), vcij(3),gij(3),vci(3),vcj(3);

muij = mi*mj/(mi+mj);

for (int i = 0; i < 3; i++) {
  vci[i] = vi[i] + 0.5*dt/mi*fi[i]*KCALMOLANGAMU_TO_ANGFS2; // v'1(dt)
  vcj[i] = vj[i] + 0.5*dt/mj*fj[i]*KCALMOLANGAMU_TO_ANGFS2; // v'2(dt)
  rij[i] = ri[i] - rj[i];                                   // r12(dt)
  vcij[i] = vci[i] - vcj[i];                                // v'12(dt)
} 
double rijvcij = rij[0]*vcij[0] + rij[1]*vcij[1] + rij[2]*vcij[2];
double r2ij = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
double lambda;

lambda = muij*rijvcij/r2ij;

// update velocities after apply second constrained
for (int i = 0; i < 3; i++) {
  vi[i] = vci[i] - lambda/mi*rij[i]; // v1(dt)
  vj[i] = vcj[i] + lambda/mj*rij[i]; // v2(dt)	
}

}
