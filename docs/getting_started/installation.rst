Installation
============

Download the `MassCCS <https://github.com/cces-cepid/massccs>`_ or clone the repository on your computer::

    git clone https://github.com/cces-cepid/massccs.git


Required Software
~~~~~~~~~~~~~~~~~

MassCCS depends on the following software:

* C++ 9.3+
* CMake 3.13+

On Ubuntu/Debian, you can simply run the following commands to install them with the package manager::

  sudo apt install gcc
  sudo apt-get install cmake

Installing
~~~~~~~~~~
.. _label:

This software uses CMake 3.13+ as its build system generator. On your terminal, run the following commands to compile the source code::
    
    cd massccs
    mkdir build  # Create build directory    
    cd build    
    cmake ..     # Generate Makefiles
    make         # Compile the program

There are additional options than can be used to compile the code. To see them, just run cmake ``cmake -LH ..``
