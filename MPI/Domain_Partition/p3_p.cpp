/*HPC-1, Fall-2016, Assignment-3, Siddhant S. Aphale, Person # - 50164327
 MPI Program to estimate the area of Mandelbrot Set
 MPI_COMM_WORLD partitioning in two groups, 
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

  //Set MPI Status
  MPI_Status status;

  //Initializing the step size in the X and Y direction in COmplex Plane
  double hx, hy;
  double xmax = 2.0, xmin = -2.0;
  double ymax = 2.0, ymin = -2.0;

  //Initializing the real and imaginary parts of Complex Numbers z and C
   
  double* x;
  double* y;

  x = new double[N];
  y = new double[N/2];
  //C = new double[N*N];

  double global_area;

  double** Creal;
  double** Cimag;

  //Splitting the rows and columns of the complex grid
  int rows = N/2;
  int columns = N;
  
  Creal = new double*[rows];
  Cimag = new double*[N/2];

  for (int i = 0; i < rows; i++)
    {
      Creal[i] = new double[columns];
      Cimag[i] = new double[columns];
    }

  double t1, t2;

  hx = (xmax - xmin)/(N - 1);
  hy = (ymax - ymin)/(N - 1);
  
  int nproc, my_proc, ierr;
  int global_id, color, new_id, new_nodes;
  MPI_Comm New_Comm;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
  color = my_proc%2;
  global_id = my_proc;
  int tag = 999;
  MPI_Comm_split(MPI_COMM_WORLD, color, my_proc, &New_Comm);
  MPI_Comm_rank(New_Comm, &new_id);
  MPI_Comm_size(New_Comm, &new_nodes);
  
  t1 = MPI_Wtime();
  if (color == 0)
    {
      double total_sum_p1 = 0.0, total_sum_p2 = 0.0, sum_p1 = 0.0, sum_p2 = 0.0;
      
      //Generating the grid points in x and y direction on the complex plane
      x[0] = -2.0;
      y[0] = -2.0;
      
      for (int i = 1; i < N; i++)
	{
	  x[i] = x[i-1] + hx;
	}
      
      for (int i = 1; i < N/2; i++)
	{
	  y[i] = y[i-1] + hy;
	}
      
      //Meshing the grid
      for (int i = 0; i < rows; i++)
	{
	  for (int j = 0; j < columns; j++)
	    {
	      Creal[i][j] = x[j];
	      Cimag[i][j] = y[i];
	    }
	}
      
      int start, end;
      
      start = ((N*(N/2))/new_nodes)*new_id;
      end = ((N*(N/2))/new_nodes)*(new_id + 1);
      
      /*for (int i = 0; i < N; i++)
	{
	for (int j = 0; j < N; j++)
	{
	C[i*N+j] = Creal[i][j];
	}
	}*/
      
      
      double r, im;
      
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
		  sum_p1 += hx*hy;
		}
	      z0 = z1;
	    }
	}
      
      ierr = MPI_Reduce(&sum_p1, &total_sum_p1, 1, MPI_DOUBLE, MPI_SUM, 0, New_Comm);
      
      if (new_id == 0)
	{
	  std::cout << "Sum1 of individual processor = " << total_sum_p1 << std::endl;
	}
      
      if (new_id == 0 && my_proc == 0)
	{
	  MPI_Recv(&total_sum_p2, 1, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &status);
	  
	  global_area = total_sum_p1 + total_sum_p2;
	  std::cout << "Global sum = " << global_area << std::endl;
	}     
      
      delete[] x;
      delete[] y;
      
      
      for (int i = 0; i < rows; i++)
	{
	  delete[] Creal[i];
	  delete[] Cimag[i];
	}
      delete[] Creal;
      delete[] Cimag;
    }
  else
    {
      double total_sum_p1 = 0.0, total_sum_p2 = 0.0, sum_p1 = 0.0, sum_p2 = 0.0;
      
      //Generating the grid points in x and y direction on the complex plane
      x[0] = -2.0;
      y[(N/2)-1] = 2;
      
      for (int i = 1; i < N; i++)
	{
	  x[i] = x[i-1] + hx;
	}
      for (int i = ((N/2)-2); i > -1; i--)
	{
	  y[i] = y[i+1] - hy;
	}
      
      //Meshing the grid
      for (int i = 0; i < rows; i++)
	{
	  for (int j = 0; j < columns; j++)
	    {
	      Creal[i][j] = x[j];
	      Cimag[i][j] = y[i];
	    }
	}
      
      int start, end;
      
      start = ((N*(N/2))/new_nodes)*new_id;
      end = ((N*(N/2))/new_nodes)*(new_id + 1);
      
      double r, im;
      
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
		  sum_p2 += hx*hy;
		}
	      z0 = z1;
	    }
	}
      
      ierr = MPI_Reduce(&sum_p2, &total_sum_p2, 1, MPI_DOUBLE, MPI_SUM, 0, New_Comm);

      if (new_id == 0 && my_proc == 1)
	{
	  MPI_Send(&total_sum_p2, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
      
      delete[] x;
      delete[] y;
      
      for (int i = 0; i < rows; i++)
	{
	  delete[] Creal[i];
	  delete[] Cimag[i];
	}
      delete[] Creal;
      delete[] Cimag;
    }
  
  t2 = MPI_Wtime() - t1;
  
  MPI_Finalize();
  
  if (my_proc == 0)
    {
      std::string file = "nproc";
      std::ofstream area_output;
      std::stringstream proc1;
      proc1 << nproc;
      
      std::string filename = "mandelbrotpart"+proc1.str()+".dat";
      area_output.open(filename.c_str());
      area_output << nproc << "\t" << t2 << "\t" << global_area << "\n";
      area_output.close();
      proc1.str("");
    }
  
  return 0;
}  
