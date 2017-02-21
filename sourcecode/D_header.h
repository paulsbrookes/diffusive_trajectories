#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>

#define Power pow 
#define Conjugate conj

using namespace std;

void D1 (complex<double>* psinew, complex<double>* psi, void *p);
void D2(complex<double>* psinew, complex<double>* psi, void *p);
double Sqrt (double x);
complex<double> Complex( double re, double im);
