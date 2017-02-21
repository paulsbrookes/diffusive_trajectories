#include "D_header.h"
#include "parameters_struct.h"
#include <complex>

using namespace std;

void  D1(complex<double>* psinew, complex<double>* psi, void *p)
{
    struct myfunc_params *params = (struct myfunc_params*)p;
    double kappa = params->kappa;
    double gamma = params->gamma;
    double drive = params->drive;
    double dq = params->dq;
    double dc = params->dc;
    double g = params->g;
    double chi = params->chi;
    int size   = params->size;
    double root2 = sqrt(2);
    double rootg = sqrt(gamma);
    double rootk = sqrt(kappa);
    #include "expectationk.dat"
    #include "expectationg.dat"
    complex<double> square_expectationk = Power(expectationk,2);
    complex<double> square_expectationg = Power(expectationg,2);
    #include "D1.dat"
}

void  D2_kappa(complex<double>* psinew, complex<double>* psi, void *p)
{
    struct myfunc_params *params = (struct myfunc_params*)p;
    double kappa = params->kappa;
    double drive = params->drive;
    double dq = params->dq;
    double dc = params->dc;
    double g = params->g;
    int size = params->size;
    double root2 = sqrt(2);
    double rootk = sqrt(kappa);
    #include "expectationk.dat"
    #include "D2_kappa.dat" 
}

void  D2_gamma(complex<double>* psinew, complex<double>* psi, void *p)
{
    struct myfunc_params *params = (struct myfunc_params*)p;
    double gamma = params->gamma;
    double drive = params->drive;
    double dq = params->dq;
    double dc = params->dc;
    double g = params->g;
    int size = params->size;
    double root2 = sqrt(2);
    double rootg = sqrt(gamma);
    #include "expectationg.dat"
    #include "D2_gamma.dat" 
}
