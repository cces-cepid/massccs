.. massccs documentation documentation master file, created by
   sphinx-quickstart on Tue Feb 21 21:30:23 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 4
   :caption: Content
   :hidden:

   self
   getting_started   
   useful

MassCCS
=======

Massive Collision Cross Section calculations for large macromolecular assemblies (MassCCS). MassCCS is an efficient parallel software
for calculating Collision Cross Section (CCS) under the trajectory method (TM) for any macromolecular structure regardless of its size, 
shape and surface rugosity, without significant loss of accuracy or performance. 

MassCCS relies on the linked-cell algorithm for computing intermolecular forces, ellipsoid projection approximation and efficient
parallelization techniques that greatly enhance its performance in comparison to available TM CCS codes.

MassCCS is written in C++ and supports OpenMP. Extensive tests on calculation accuracy, speed up gains, and scalability with
system size were performed. MassCCS efficiency is particularly scalable for multiple CPU cores.

If you use MassCCS in your research please, cite the following papers:

* `S. Cajahuaringa, D. L. Z. Caetano, L. N. Zanotto, G. Araujo and M. S. Skaf, MassCCS: A high performance collision cross-section software for large macromolecules, J. Chem. Inf. Model. 2023, 63, 11, 3557–3566. <https://doi.org/10.1021/acs.jcim.3c00405>`_

* `S. Cajahuaringa, L. N. Zanotto, D. L. Z. Caetano, S. Rigo, H. Yviquel, M. S. Skaf and G. Araujo, Ion-Molecule Collision Cross-Section Simulation using Linked-cell and Trajectory Parallelization, 2022 IEEE 34th International Symposium on Computer Architecture and High Performance Computing (SBAC-PAD), Bordeaux, France, 2022, pp. 150-159. <https://ieeexplore.ieee.org/document/9980906>`_

============
Availability
============

MassCCS can be :doc:`easily installed </getting_started/installation>`  using package manager.

Source code is freely available at `GitHub <https://github.com/cces-cepid/massccs>`_ under the LGPL-3.0 license.

========
Contacts
========

Please `submit a github issue <https://github.com/cces-cepid/massccs>`_ to report bugs and to request new features.

Alternatively you may find the developer `here <mailto:samuelcm@unicamp.br>`_ . Please visit our `website <https://github.com/cces-cepid/massccs>`_ for more information.
