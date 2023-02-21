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

#ifndef MASSCCS_V0_CONSTANTS_H
#define MASSCCS_V0_CONSTANTS_H

#include <cmath>

#define EMPTY (-1)

// atomic mass unit [kg]
#define AMU_TO_KG (1.66053904E-27)

// Mol
#define MOL (6.02214078E23)

// Boltzmann Constant [J/K]
#define BOLTZMANN_K (1.38064E-23)

// Calorie to Joule
#define CAL_TO_J (4.184)

// Angstrom to m
#define ANG_TO_M (1.0E-10)

// Femtosecond to second
#define FS_TO_S (1.0E-15)

// degree to radians
#define DEGREE_TO_RAD (M_PI / 180.0)

// radians to degree
#define RAD_TO_DEGREE (180.0 / M_PI)

// Force from kcal/mol/Ang/amu to Ang/fs^2
#define KCALMOLANGAMU_TO_ANGFS2 (1.0/48.88821291/48.88821291)
//#define KCALMOLANGAMU_TO_ANGFS2 4.1840000545735656E-4
// Joule to eV
#define J_TO_eV (6.242E18)
// ev to kcal/mol
#define eV_TO_KCAL_MOL (23.061)

// alpha to kcal/mol
#define ALPHA_TO_KCAL_MOL (331.842254885)

// omega in fs⁻¹
#define OMEGA_TO_FS_INV (2.45400505E8)

// coulomb constant in kcal/mol
#define KCOUL (332.06348078020324)

// amu.Ang2/fs2 to kcal/mol
#define amuAngfs2_to_KCAL_MOL (2390.07)

// Input default parameters
#define NPROBE 10000 
#define NITER 10
#define SEED 20162104   
#define TIMESTEP 10.0
#define TEMPERATURE 298.0
#define SKIN 0.01
#define SHORT_CUTOFF 12.0
#define INNER_SHORT_CUTOFF 10.0
#define LONG_CUTOFF 25.0
#define INNER_LONG_CUTOFF 22.0
#define ALPHA_HE 0.204956  
#define ALPHA_N2 1.710
#define ALPHA_CO2 2.911

#endif // MASSCCS_V1_CONSTANTS_H
