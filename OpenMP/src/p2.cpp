//HPC-1 Assignment-2, Siddhant Aphale, Person # - 50164327
//Problem 2 -  Mandelbrot Set Generation and Area Calculation - Serial Code

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <vector>
#include "ComplexNumber.hpp"
#include "timer.h"
#include <papi.h>

using namespace std;

int main(int argc, char* argv[])
{
  //Number of nodes
  vector<int> N{100, 500, 1000, 1500};

  //Initializing the step size in the X and Y direction in COmplex Plane
  double hx, hy;
  double xmax = 2.0, xmin = -2.0;
  double ymax = 2.0, ymin = -2.0;
  timespec before, after, time_diff;


  for (int l = 0; l < N.size(); l++)
    {
      //Initializing the real and imaginary parts of Complex Numbers z and C
      vector<double> x(N[l], 0.0);
      vector<double> y(N[l], 0.0);
      vector<double> C(N[l]*N[l] , 0.0);
      vector<double> time_s(N[l], 0.0);

      vector<vector<double> > Creal(N[l], vector<double>(N[l], 0.0));
      vector<vector<double> > Cimag(N[l], vector<double>(N[l], 0.0));

      /*float real_time, proc_time, mflops;
      long long flpops;
      float ireal_time, iproc_time, imflops;
      long long iflpops;*/
      int retval, w;
      long long values[2];
      int PAPI_events[] = {PAPI_L1_TCM, PAPI_L2_TCM};

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

      double r = 0.0, im = 0.0, s = 0.0;

      /*if ((retval = PAPI_flops(&ireal_time, &iproc_time, &iflpops, &imflops)) < PAPI_OK)
	{
	  std::cout << "Could not initialize PAPI_flops \n";
	  std::cout << "retval : " << retval << "\n";
	  return(1);
	}*/
      
      get_time(&before);

      PAPI_library_init(PAPI_VER_CURRENT);
      w = PAPI_start_counters(PAPI_events, 2);
      for (int i = 0; i < N[l]*N[l]; i++)
	{
	  int p = i/N[l];
	  int q = i%N[l];
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
		  s += hx*hy;
		}
	      z0 = z1;
	    }
	}
      get_time(&after);

      /*if ((retval = PAPI_flops(&real_time, &proc_time, &flpops, &mflops))<PAPI_OK)
	{
	  std::cout << "retval : " << retval << "\n";
	  return(1);
	}*/

      PAPI_read_counters(values, 2);

      diff(&before,&after,&time_diff);

      // Time in seconds
      time_s[l] = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;


      cout << "Time required = " << time_s[l] << "\n";
      cout << "Area of the Mandelbrot set = " << s << "\n";

      //std::cout << "Real time : " << real_time << "\tProc_time : " << proc_time << "\tTotal flpops : " << flpops << "\tMFLOPS : " << mflops << "\n";

      std::cout << "L1 cache misses = " << values[0] << "L2 cache misses = " << values[1] << "\n";

    }

  return 0;
}  
