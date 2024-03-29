# MassCSS data set

This directory contains:

* The files for all proteins considered in this work are founded in the following ["link"](https://drive.google.com/file/d/1NYWqjpee3eDW_63tvzY9bVBhLHLGi3jX/view?usp=share_link) (in extended XYZ-format), for more information see the table S2 of supplementary material

* A template of input.json file using in this work, it is only necessary to change the buffer gas type (He or N2) and the target filename (protein file), for more information review the ["MassCCS doc"](https://massccs.readthedocs.io/en/latest/).

* The outputs directory contains the ouput files for the CCS calculation for Helium and Nitorgen buffer gases, run in a single node with 1, 2, 4, 8 and 16 threads.

* The file ccs_of_pdb.csv contains CCS values calculated by MassCCS and IMPACT for over 110,000 structures drawn from the ["Protein Data Bank"](https://www.rcsb.org), the structure of the table is (sorted by weigth mass): 

| index | id | file |	natoms | mass | charge | ccs_He | err_He | ccs_N2 | err_N2 | ccs_PA | ratio_He | ratio_N2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1165 | 1ay3 | pdb1ay3.pqr |	25 | 0.1731 | -7.77156E-16 | 90.4027 | 4.26514 | 141.957 | 2.5643 |	82.9363 | 1.09002 | 1.71163 |

where: **id** is the PBD id of the protein, **file** is the associated pqr file of the protein generated by the ["PDB2PQR"](https://github.com/Electrostatics/pdb2pqr) software at pH = 7.0, **natoms** is the number of atoms for the protonated protein, **mass** is the weigth mass in kDa, **charge** is liquid charge in e, **ccs_He** and **ccs_N2** are the CCS (in A²) for Helium and Nitrogen buffer gases using MassCCS, **err_He** and **err_N2** are the associated errors, **ccs_PA** is the CCS using ["IMPACT"](https://process.innovation.ox.ac.uk/software/p/10126c/(impact---commercial-use-only)/1) software, **ratio_He** and **ratio_N2** are the ratio between the CCS using the MassCCS and IMPACT.
