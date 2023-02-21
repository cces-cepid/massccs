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

#include "headers/GasBuffer.h"

GasBuffer::GasBuffer(int gas_buffer_flag) {

if (gas_buffer_flag == 1) {
  natoms = 1;
  gas_type = "He";
} else if (gas_buffer_flag == 2) {
  natoms = 3; 
  gas_type = "N2";
} else if (gas_buffer_flag == 3) {
  natoms = 3; 
  gas_type = "CO2";
}

x = new double[natoms];
y = new double[natoms];
z = new double[natoms];
vx = new double[natoms];
vy = new double[natoms];
vz = new double[natoms];
q = new double[natoms];
m = new double[natoms];
eps = new double[natoms];
sig = new double[natoms];
atomName = new string[natoms];

gas_properties();

rcm[0] = 0.0;
rcm[1] = 0.0;
rcm[2] = 0.0;

}

void GasBuffer::gas_properties() {

if (gas_type == "He") {
  x[0] = 0.0;
  y[0] = 0.0;
  z[0] = 0.0;
  vx[0] = 0.0;
  vy[0] = 0.0;
  vz[0] = 0.0;
  q[0] = 0.0;
  m[0] = 4.0026;
  eps[0] = 0.0309008;
  sig[0] = 3.043;
  atomName[0] = "He";
  mass = m[0];
  d = 0.0;
} else if (gas_type == "N2") {
  x[0] = 0.0;
  y[0] = 0.0;
  z[0] = 0.5488;
  x[1] = 0.0;
  y[1] = 0.0;
  z[1] = -0.5488; 
  x[2] = 0.0;
  y[2] = 0.0;
  z[2] = 0.0;
  vx[0] = 0.0;
  vy[0] = 0.0;
  vz[0] = 0.0;
  vx[1] = 0.0;
  vy[1] = 0.0;
  vz[1] = 0.0;
  vx[2] = 0.0;
  vy[2] = 0.0;
  vz[2] = 0.0;
  q[0] = -0.4825;
  q[1] = -0.4825;
  q[2] = 0.965; 
  m[0] = 14.007;
  m[1] = 14.007;
  m[2] = 0.0;
  eps[0] = 1.0;
  eps[1] = 1.0;
  eps[2] = 0.0;
  sig[0] = 0.0;
  sig[1] = 0.0;
  sig[2] = 0.0;
  atomName[0] = "N_1";
  atomName[1] = "N_2";
  atomName[2] = "Dummy";
  mass = m[0] + m[1];
  d = abs(z[0]-z[1]);
} else if (gas_type == "CO2") {
  // force field paramter of Harris and Yung (1995)
  // oxygen 1
  x[0] = 0.0;
  y[0] = 0.0;
  z[0] = -1.49;
  // carbon
  x[1] = 0.0;
  y[1] = 0.0;
  z[1] = 0.0;
  // oxygen 2
  x[2] = 0.0;
  y[2] = 0.0;
  z[2] = 1.49;
  vx[0] = 0.0;
  vy[0] = 0.0;
  vz[0] = 0.0;
  vx[1] = 0.0;
  vy[1] = 0.0;
  vz[1] = 0.0;
  vx[2] = 0.0;
  vy[2] = 0.0;
  vz[2] = 0.0;
  q[0] = -0.3256;
  q[1] = 0.6512;
  q[2] = -0.3256;
  m[0] = 15.999;
  m[1] = 12.011;
  m[2] = 15.999;
  eps[0] = 0.159; 
  eps[1] = 0.055; 
  eps[2] = 0.159; 
  sig[0] = 3.033; 
  sig[1] = 2.757; 
  sig[2] = 3.033; 
  atomName[0] = "O1";
  atomName[1] = "C1";
  atomName[2] = "O1";
  mass = m[0] + m[1] + m[2];
  d = abs(z[0]-z[2]);
} 

} 
