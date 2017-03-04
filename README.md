# Mandelbrot_Set_OpenMP_MPI
C++ Implementation of Mandelbrot Set Area Calculation with OpenMP and MPI


OpenMP:
p2.cpp is the serial code with PAPI Profiling
p2dynamic.cpp, p2guided.cpp, p2static.cpp are the respective scheduling stategies of OpenMP
ComplexNumber.hpp is the class that generates Complex Numbers, calculates modulus and arguements

MPI:
Domain was partitioned by creating a new communicator and is under "Domain_Partition"
The code without domain partition is under "No_Partition"
