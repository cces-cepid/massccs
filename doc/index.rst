.. MassCCS documentation documentation master file, created by
   sphinx-quickstart on Tue Oct 11 13:02:18 2022.
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

* :ref:`S. Cajahuaringa, D. L. Z. Caetano, L. N. Zanotto, G. Araujo and M. S. Skaf, MassCCS: A high performance collision cross-section software for large macromolecules, submitted to Analytical Chemistry`

* :ref:`S. Cajahuaringa, L. N. Zanotto, D. L. Z. Caetano, S. Rigo, H. Yviquel, M. S. Skaf and G. Araujo, Ion-Molecule Collision Cross-Section Simulation using Linked-cell and Trajectory Parallelization, IEEE SBAC-PAD 2022`

============
Availability
============

MassCCS can be :doc:`easily installed </getting_started/installation>`  using package manager.

Source code is freely available at `GitHub <https://github.com/cepid-cces/massccs>`_ under the LGPL-3.0 license.

========
Contacts
========

Please `submit a github issue  <https://github.com/cepid-cces/massccs>`_ to report bugs and to request new features.

Alternatively you may find the developer `here <mailto:skaf@unicamp.br>`_ . Please visit our `website <https://github.com/cepid-cces/massccs>`_ for more information.
