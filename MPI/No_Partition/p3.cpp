/*HPC-1, Fall-2016, Assignment-3, Siddhant S. Aphale, Person # - 50164327
 MPI Program to estimate the area of Mandelbrot Set
 Study of Strong and Weak Scaling*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include "ComplexNumber.hpp"


int main(int argc, char* argv[])
{
  //Number of nodes
  int N = 1200;
  //int z = sizeof(N)/sizeof(N[0]);

  //Initializing the step size in the X and Y direction in COmplex Plane
  double hx, hy;
  double xmax = 2.0, xmin = -2.0;
  double ymax = 2.0, ymin = -2.0;


  //for (int l = 0; l < N; l++)
  //{
      //Initializing the real and imaginary parts of Complex Numbers z and C

  double* x;
  double* y;
  double* C;

  x = new double[N];
  y = new double[N];
  C = new double[N*N];

  double** Creal;
  double** Cimag;
  Creal = new double*[N];
  Cimag = new double*[N];

  for (int i = 0; i < N; i++)
    {
      Creal[i] = new double[N];
      Cimag[i] = new double[N];
    }

      //Generating the grid points in x and y direction on the complex plane
      x[0] = -2.0;
      y[0] = -2.0;

      hy = (ymax - ymin)/(N - 1);

      for (int i = 1; i < N; i++)
	{
	  x[i] = x[i-1] + hy;
	  y[i] = y[i-1] + hy;
	}

      //Meshing the grid
      for (int i = 0; i < N; i++)
	{
	  for (int j = 0; j < N; j++)
	    {
	      Creal[i][j] = x[j];
	      Cimag[i][j] = y[i];
	    }
	}

      for (int i = 0; i < N; i++)
	{
	  for (int j = 0; j < N; j++)
	    {
	      C[i*N+j] = Creal[i][j];
	    }
	}

      unsigned int a;
      double r = 0.0, im = 0.0;
      long int start, end;
      int nproc, my_proc, ierr;

      double sum, total_sum, area;
      double t1, t2;

      MPI_Init(&argc, &argv);
      MPI_Comm_size(MPI_COMM_WORLD, &nproc);
      MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);

      if (my_proc == 0)
	{
	  std::cout << "MPI tasks = " << nproc << "\n";
	  std::cout << "Number of node points = " << N << "\n";
	}

      if (a = (N*N%nproc))
	{
	  if (my_proc == 0)
	    {
	      std::cout << "Number of processors must divide exactly" << N << std::endl;
	      MPI_Finalize();
	      return(1);
	    }
	}

      t1 = MPI_Wtime();
      start = ((N*N)/nproc)*my_proc;
      end = ((N*N)/nproc)*(my_proc+1);
      sum = 0.0;
      total_sum = 0.0;

      hx = (xmax - xmin)/(N - 1);
      
    
      for (int i = start; i < end; i++)
	{
	  int p = i/N;
	  int q = i%N;
	  r = Creal[p][q];
	  im = Cimag[p][q];

	  ComplexNumber z1;
	  ComplexNumber c0(r, im);
	  ComplexNumber z0;

	  for (int j = 0; j < 10000; j++)
	    {
	      z1 = z0.Power(2) + c0;

	      if (z1.Modulus() > 2)
		{
		  break;
		}
	      else if ( j == 9999)
		{
		  sum += hx*hy;
		}
	      z0 = z1;
	    }
	}

      ierr = MPI_Reduce(&sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      t2 = MPI_Wtime() - t1;

      area = total_sum;

      if (my_proc == 0)
	{
	  std::cout << "Area of Manelbrot Set = " << area << "\n";
	  std::cout << "Time required = " << t2 << "\n";
	}

      MPI_Finalize();

      if (my_proc == 0)
	{
	  std::string file = "nproc";
	  std::ofstream area_output;
	  std::stringstream proc1;
	  proc1 << nproc;

	  std::string filename = "mandelbrot"+proc1.str()+".dat";
	  area_output.open(filename.c_str());
	  area_output << nproc << "\t" << t2 << "\t" << area << "\n";
	  area_output.close();
	  proc1.str("");
	}

      delete[] x;
      delete[] y;
      delete[] C;

      for (int i = 0; i < N; i++)
	{
	  delete[] Creal[i];
	  delete[] Cimag[i];
	}
      delete[] Creal;
      delete[] Cimag;
      


  return 0;
}  
