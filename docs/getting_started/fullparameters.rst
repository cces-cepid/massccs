Complete Parameters description 
===============================

Each possible program option that can be used in the parameter input file is listed below along with a brief description of its functionality. It is worth mentioning that not all CHAPTER options need to be included in the parameter input file. The user can choose which ones (if any) are appropriate for their simulation. In this way, all options, including the parameter input file itself, are optional. As a result, each changeable variable also has a default value.

numberProbe-tag
---------------
Number of trajectories calculations per CCS integral.

Default value: 1000

nIter-tag
---------
Number of CCS integrals.

Default value: 10

nthreads-tag
------------
Number of threads allocated for execution.

Default value: maximal available threads

Temp-tag
--------
Temperature, in Kelvin, of the buffer gas.

Default value: 298.0 K.

seed-tag
--------
The seed of random number generator is used to positon the atoms/molecules of the buffer gas around the target molecule and generate their velocity distribution according to the chosen temperature

dt-tag
------
Time step, in femtosecods (fs).

Default value: 10 fs

skin-tag
--------
Skin size of unit cell used in the Linked-Cell algorithm.

Default value: 0.01 Å

GasBuffer-tag
-------------
Buffer gas selection - helium (He), nitrogen (N2), or carbon dioxide (CO2).

Default value: He.

Equipotential-tag
-----------------
Ellipsoid type selection. If selected "yes", an equipotential surface is calculated according to the shape of the target molecule. Then, an ellipsoid is fitted to this surface and magnified uniformly in each direction until the equipotential surface fits completely inside the ellipsoid. This procedure is quite time consuming for large systems. If "no" is selected, the calculation of the ellipsoid is based solely on the geometric shape of the target molecule. 

Defaltu value: no. 

Short-range cutoff-tag
----------------------
If sected "yes", a cutoff radius is applied when calculating Lennard-Jones forces. If "no" is selected, the calculation of Lennard-Jones forces is performed on all atoms. 

Default value: yes

LJ-cutoff-tag
-------------

Cutoff radius of Lennard-Jones forces.

Default value: 12.0 Å

Long-range forces-tag
---------------------
If selected "yes", long-range forces (ion-quadrupole forces) are also calculated. If "no" is selected, only Lennard-Jones (and ion-induced dipole) forces are calculated.

Default value: no


Long-range cutoff-tag
---------------------
If sected "yes", a cutoff radius is applied when calculating long-range forces. If "no" is selected, the calculation of long-range forces is performed on all atoms.

Default value: no

Coul-cutoff-tag
---------------
Cutoff radius of long-range forces.

Default value: 25.0 Å

polarizabilty-tag
-----------------
If selected "yes", ion-induced dipole forces are also calculated. If "no" is selected,
only Lennard-Jones (and ion-quadrupole) forces are calculated.

Default value: no

.. warning:: 
    MassCCS can also handle other linear buffer molecules like CO2, although the 
    interaction parameters have not been parameterized for CCS calculations using carbon dioxide.

