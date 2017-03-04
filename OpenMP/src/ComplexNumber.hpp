//Header File to generate complex numbers and perform various operations
//Calculate Modulus, Argument and Power of Complex Number

#ifndef COMPLEXNUMBERHEADERDEF
#define COMPLEXNUMBERHEADERDEF

#include <iostream>
#include <cmath>

class ComplexNumber
{
private:
  double mReal;
  double mImaginary;
public:
  ComplexNumber()
  {
    mReal = 0.0;
    mImaginary = 0.0;
  }

  ComplexNumber(double x , double y)
  {
    mReal = x;
    mImaginary = y;
  }

  double Modulus() const
  {
    return sqrt(mReal*mReal + mImaginary*mImaginary);
  }

  double Argument() const
  {
    return atan2(mImaginary, mReal);
  }

  ComplexNumber Power(double n) const
  {
    double mod = Modulus();
    double arg = Argument();
    double mod_result = pow(mod, n);
    double arg_result = arg*n;
    double real = mod_result*cos(arg_result);
    double imaginary = mod_result*sin(arg_result);
    ComplexNumber z(real, imaginary);
    return z;
  }

  ComplexNumber operator = (const ComplexNumber z)
  {
    mReal = z.mReal;
    mImaginary = z.mImaginary;
    return *this;
  }

  ComplexNumber operator + (const ComplexNumber& z) const
  {
    ComplexNumber w;
    w.mReal = mReal + z.mReal;
    w.mImaginary = mImaginary + z.mImaginary;
    return w;
  }

};

#endif
