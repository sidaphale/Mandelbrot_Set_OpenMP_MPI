//HPC-1 Assignment-2, Siddhant Aphale, Person # - 50164327
//Problem 2 -  Mandelbrot Set Generation and Area Calculation - OpenMP dynamic

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>
#include "ComplexNumber.hpp"
#include "timer.h"
#include <omp.h>

int main(int argc, char* argv[])
{
  //Number of nodes
  std::vector<int> N{500};

  //Initializing the step size in the X and Y direction in COmplex Plane
  double hx, hy;
  double xmax = 2.0, xmin = -2.0;
  double ymax = 2.0, ymin = -2.0;
  timespec before, after, time_diff;


  for (int l = 0; l < N.size(); l++)
    {
      //Initializing the real and imaginary parts of Complex Numbers z and C
      std::vector<double> x(N[l], 0.0);
      std::vector<double> y(N[l], 0.0);
      std::vector<double> C(N[l]*N[l] , 0.0);
      std::vector<double> time_s(16, 0.0);

      std::vector<std::vector<double>> Creal(N[l], std::vector<double>(N[l], 0.0));
      std::vector<std::vector<double>> Cimag(N[l], std::vector<double>(N[l], 0.0));

      //Defining the step size
      hx = (xmax-xmin)/(N[l]-1);
      hy = (ymax-ymin)/(N[l]-1);

      //Generating the grid points in x and y direction on the complex plane
      x[0] = -2.0;
      y[0] = -2.0;

      for (int i = 1; i < N[l]; i++)
	{
	  x[i] = x[i-1] + hx;
	  y[i] = y[i-1] + hy;
	}

      //Meshing the grid
      for (int i = 0; i < N[l]; i++)
	{
	  for (int j = 0; j < N[l]; j++)
	    {
	      Creal[i][j] = x[j];
	      Cimag[i][j] = y[i];
	    }
	}

      for (int i = 0; i < N[l]; i++)
	{
	  for (int j = 0; j < N[l]; j++)
	    {
	      C[i*N[l]+j] = Creal[i][j];
	    }
	}

      for (int t=0; t < 16; t++)
	{
	  

	  double r = 0.0, im = 0.0, s = 0.0;
	  int i, p, q,j;
	  omp_set_num_threads(t);
	  get_time(&before);
	  
#pragma omp parallel for private(i,p,q,j) reduction(+:s) schedule(dynamic, 50)
	  for (i = 0; i < N[l]*N[l]; i++)
	    {
	      p = i/N[l];
	      q = i%N[l];
	      r = Creal[p][q];
	      im = Cimag[p][q];
	      
	      ComplexNumber z1;
	      ComplexNumber c0(r, im);
	      ComplexNumber z0;
	      
	      for (j = 0; j < 10000; j++)
		{
		  z1 = z0.Power(2) + c0;
		  
		  if (z1.Modulus() > 2)
		    {
		      break;
		    }
		  else if ( j == 9999)
		    {
		      s += hx*hy;
		    }
		  z0 = z1;
		}
	    }
	  get_time(&after);
	  
	  diff(&before,&after,&time_diff);
	  
	  // Time in seconds
	  time_s[t] = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;
	  std::cout << "Area of Mandelbrot set = " << s << "\n";
	}

      std::ofstream write_dynamic("dynamic.dat");
      write_dynamic.precision(16);
      for (int i = 0; i < 16; i++)
	{
	  write_dynamic << i << "\t" << time_s[i] << "\n";
	}      
    }
  
  return 0;
}  
