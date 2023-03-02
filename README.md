# MassCSS: A high performance collision cross-section software

Massive Collision Cross Section calculations for large macromolecular assemblies (MassCCS). MassCCS is an efficient parallel software for calculating Collision Cross Section (CCS) under the trajectory method (TM) for any macromolecular structure regardless of its size, shape and surface rugosity, without significant loss of accuracy or performance.

MassCCS relies on the linked-cell algorithm for computing intermolecular forces, ellipsoid projection approximation and efficient parallelization techniques that greatly enhance its performance in comparison to available TM CCS codes.

MassCCS is written in C++ and supports OpenMP. Extensive tests on calculation accuracy, speed up gains, and scalability with system size were performed. MassCCS efficiency is particularly scalable for multiple CPU cores.

If you use MassCCS in your research please, cite the following papers:

["S. Cajahuaringa, D. L. Z. Caetano, L. N. Zanotto, G. Araujo and M. S. Skaf, MassCCS: A high performance collision cross-section software for large macromolecular assemblies, in preparation"]()

["S. Cajahuaringa, L. N. Zanotto, D. L. Z. Caetano, S. Rigo, H. Yviquel, M. S. Skaf and G. Araujo, Ion-Molecule Collision Cross-Section Simulation using Linked-cell and Trajectory Parallelization, 2022 IEEE 34th International Symposium on Computer Architecture and High Performance Computing (SBAC-PAD), Bordeaux, France, 2022, pp. 150-159."](https://ieeexplore.ieee.org/abstract/document/9980906)

## Install

This software uses CMake 3.13+ as its build system generator. On your terminal,
run the following commands to compile the source code:

```bash
cd massccs
mkdir build # Create build directory
cd build
cmake .. # Generate Makefiles
make  # Compile the program
```

There are optional compilation option that can be used, to see a list of all of them run `cmake -LH ..`.

## Run

```bash
cd .. # need to execute massccs from root dir
./build/massccs ./input.json
```

For more information about the software installation and used review the ["MassCCS doc"](https://massccs-doc.readthedocs.io/en/latest/)

Author & Contact:
--------------
Samuel Cajahuaringa - samuelcajahuaringa@gmail.com

