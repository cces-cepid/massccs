CCS calculation
===============

Structure File
~~~~~~~~~~~~~~

In this example, we calculate the CCS value for the Ubiquitin protein (PDB-ID: 1UBQ)
The structure may be downloaded from the `Protein Data Bank <https://www.rcsb.org>`_ ::

   wget https://files.rcsb.org/download/1ubq.pdb

MassCCS support molecule file in ``PQR`` and ``XYZ-like`` format. 

Before calculating the CCS by the MassCCS software, it is necessary to assign partial charges 
to the atoms of the PDB structure file. To address this issue, we can use the `PDB2PQR <http://pdb2pqr.sourceforge.net/>`_ 
software package that provides a utility for converting protein files from PDB format to PQR format. 
With no installation PDB2PQR can be used free of charge through the web server
`PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_. For more information about
PDB2PQR tool you can visit: https://pdb2pqr.readthedocs.io/en/latest/ .

For the case of ubiquitin, when we select pH=4.0 in the PDB2PQR web server and AMBER forcefield, we generate a PQR file 
where the total charge of the protein is +4e.

.. note::
    Before the CCS calculation by the MassCCS software, we recommend relax the PDB structures by means energy minimization. 

Parameter Input File
~~~~~~~~~~~~~~~~~~~~

It is necessary to use the following input parameter file, in json format, to run MassCCS: 

.. code-block:: text
   :caption: input.json

   {   
        "targetFileName": "1ubq.pqr",   # Input structure file name in PQR format
        "numberProbe" : 10000,          # number of trajectories calculations per CCS integral
        "nIter" : 10,                   # number of CCS of Monte Carlo integration
        "seed" : 2104,                  # random seed
        "dt" : 10.0,                    # time step in fs
        "Temp" : 298.0,                 # temperature of buffer gas in Kelvin
        "skin" : 0.01,                  # skin size of unit cell in Angstrom
        "GasBuffer" : "He",             # buffer gas type   
        "Equipotential" : "no",         # equipotential calculation   
        "Short-range cutoff" : "yes",   # apply cutoff on Lennard-Jones forces
        "LJ-cutoff" : 12.0,             # radius cutoff in Angstrom
        "Long-range forces" : "yes",    # apply long-range forces
        "Long-range cutoff" : "yes",    # apply cutoff on long-range forces  
        "Coul-cutoff" : 25.0,           # radius cutoff in Angstrom
        "polarizability" : "yes"        # apply polarizability
    }

Run Simulation
~~~~~~~~~~~~~~

To run the simulation execute the following command: ::

    cd .. # need to execute massccs from root dir
    ./build/massccs ./input.json

If your structure is in PQR or XYZ format is only necessary include in the input.json file and run the simulation.

Output
~~~~~~

We run this exemple (ubquitin) on CPU consisting of a Dual Intel Xeon E5-2630 v3 processor with 4 cores each and AVX2, running at 2.4 GHz.

The ouput printed on screen is:

.. code-block:: text
   :caption: simulation parameters

   *********************************************************
   INPUT:: Simulation Parameters
   *********************************************************
   target filename                  : 1ubq.pqr
   number of probe                  : 10000
   number of iterarions             : 10
   number of threads                : 4
   seed number                      : 2104
   gas buffer                       : He
   Target Temperature (K)           : 298
   timestep (fs)                    : 10
   Skin cell size (Ang)             : 0.01
   Equipotential                    : no
   Cut short-range interaction      : yes
   LJ cutoff (Ang)                  : 12
   Apply long-range interaction     : yes
   Cut long-range interaction       : yes
   Coulomb cutoff (Ang)             : 25
   Apply induced-dipole interaction : yes
   alpha (Ang^3)                    : 0.204956

.. code-block:: text
   :caption: target molecule information

   *********************************************************
   MOLECULE:: orientation around the inertia principal axis
   *********************************************************
   Inertia matrix a.u. Ang^2:
   {  812646  -103409  -125747  }
   {  -103409  803376  -61197.4  }
   {  -125747  -61197.4  758127  }
   Rotation matrix applied to molecule:
   {  0.798187  -0.433666  -0.418128  }
   {  0.59767  0.48319  0.639779  }
   {  -0.0754147  -0.760566  0.644865  }
   Molecule radius: 25.7711 Ang
   orientation time of molecule target: 0.0164136 s
   target mass: 8568.9 amu
   gas mass: 4.0026 amu
   reduce mass: 4.00073 amu
   charge state: 4 e

.. code-block:: text
   :caption: ellipsoid projection approximation

   *********************************************************
   Geometric Ellipsoid:
   *********************************************************
   maximal distances: 13.3257  15.9013  24.8827  Ang
   Initial axis length: 63.3357  65.9113  74.8927  Ang
   Ellipsoid axis length: 65.0857  67.6613  76.6427  Ang
   ellipsoid calculation time: 7.5278e-05 s
   maximal impact parameter: 76.6427 Ang

.. code-block:: text
   :caption: linked-cell list information

   *********************************************************
   Linked-cell
   ********************************************************* 
   Numbers of cells: 1716
   Nx: 11 Ny: 12 Nz: 13
   Filled cells: 29
   Empty cells: 1687
   Average atoms per cell: 42.5862
   Maximum atoms per cell: 194
   Minimum atoms per cell: 1
   Simulation box:
   lx: 132.11 Ang
   ly: 144.12 Ang
   lz: 156.13 Ang
   linked-cell calculation time: 0.0037535 s

.. code-block:: text
   :caption: trajectories information

   *********************************************************
   Trajectory calculations
   *********************************************************
   Ntraj: 10000
   Nfree: 1756
   Nscatter: 8244
   Nlost: 0
   omega: 1063.65
   Ntraj: 10000
   Nfree: 1715
   Nscatter: 8285
   Nlost: 0
   omega: 1078.6
   :
   : 
   :  
   Ntraj: 10000
   Nfree: 1738
   Nscatter: 8262
   Nlost: 0
   omega: 1135.46
   CCS time: 39.3628 s
   *********************************************************
   average value of CCS = 1126.88 Ang^2
   error value of CCS = 21.6755 Ang^2
   Total time: 39.3839 s
   Program finished...
 
