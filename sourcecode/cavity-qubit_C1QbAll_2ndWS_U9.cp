#include <chrono>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <string>
#include <unistd.h>
#include <complex>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include "D_header.h"
#include "parameters.h"
#include "parameters_struct.h"
#include <stdio.h>
#include <sstream>



gsl_rng * r; 

double photon_histogram[cavity_levels][snapshots_number];
double transmon_histogram[transmon_levels][snapshots_number];
double tnorm[bins];
int sample_cnt = 0;
int progress = 0;
int writing_in_progress = 0;
const int glb_size = cavity_levels * transmon_levels;

const double PI = 3.14159265;
const int precision = 10;

void  EvolveStep(complex<double>* psinew, complex<double>* psi, void *p)
{
   struct myfunc_params *params = (struct myfunc_params*)p;
   int size   = params->size;
   double dt  = params->dt;
   double dW  = sqrt(dt)*gsl_ran_ugaussian_ratio_method(r);



   complex<double> psitilde[size];
   complex<double> psiplus[size];
   complex<double> psiminus[size];
   complex<double> d1psi[size];
   complex<double> d2psi[size];
   complex<double> d1psitilde[size];
   complex<double> d2psiplus[size];
   complex<double> d2psiminus[size];

   D1(d1psi,psi,p);
   D2(d2psi,psi,p);

   for (int j = 0; j < size; j++){
       psitilde[j] = psi[j] + dt*d1psi[j];
       psiplus[j] = psitilde[j];
       psiminus[j] = psitilde[j];
       
       psitilde[j] += dW*d2psi[j];
       psiplus[j] += sqrt(dt)*d2psi[j];
       psiminus[j] -= sqrt(dt)*d2psi[j];
   }
   
   D1(d1psitilde,psitilde,p);
   D2(d2psiplus,psiplus,p);
   D2(d2psiminus,psiminus,p);
   
   for (int j = 0; j < size; j++){
       psinew[j] = psi[j] + 0.5*dt*(d1psitilde[j] + d1psi[j])
                   + 0.25*dW*(d2psiplus[j] + d2psiminus[j] + 2.0*d2psi[j])
                   + 0.25*(dW*dW - dt)*(d2psiplus[j] - d2psiminus[j])/sqrt(dt);
   } 
   double lnorm = 0.0;
   for (int j = 0; j < size; j++) lnorm= lnorm + pow(abs(psinew[j]),2);
   lnorm = sqrt(lnorm);
   for (int j = 0; j < size; j++) psinew[j] = psinew[j]/lnorm;
}


void TrajectoryAdd(void *p, double time, int trajectory_index, string str)
{
    struct myfunc_params *params = (struct myfunc_params*)p;

    int size   = params->size;
    double dt  = params->dt;
    
    double endtime_lcl = time;
    int bins_lcl = int(time/dt);
    
    complex <double> psi[size];
    complex <double> psinew[size];
    complex <double> trace[snapshots_number];
    complex <double> photons[snapshots_number]; //expectation value <n>
    complex <double> sigmaz[snapshots_number];
    complex <double> sigmaplus[snapshots_number];
    complex <double> sigmaminus[snapshots_number];
    complex <double> averagealpha[snapshots_number];
    complex <double> averagealphasigmaplus[snapshots_number];
    complex <double> S[snapshots_number];
    
    int snap_cnt = 0;
    complex<double> psi_snapshots[size][bins_lcl/snap_to_snap]; //declaring matrix to hold snapshots of wavefunction

    
    for (int j = 0; j < snapshots_number; j++){
        photons[j]=0.;
        sigmaz[j]=0.;
        sigmaplus[j]=0.;
        sigmaminus[j]=0.;
        averagealpha[j]=0.;
        averagealphasigmaplus[j]=0.;
	trace[j]=0.;
    }	
	
    // here comes the initialization of the state
    for ( int a = 0; a < size; a++ ){
        psi[a]=0.;
    }
    psi[transmon_levels*photon_number] = 1.0;
    psi[transmon_levels*photon_number + 1] = 0.0;

    //for ( int a = 0; a < glb_size; a++ ) {
    //    psi[a] = 0.5*sin(a + 1) + Complex(0, 1) * cos(a + 1);
    //}

    
    //set matrix elements to zero
    for( int i = 0; i < size; i++){
        for( int k = 0; k < (bins_lcl/snap_to_snap); k++ ){
            psi_snapshots[i][k] = 0.;
        }
    }   

    int mod_cnt = 0;





    int i = 1;
    for (double t = 0.0 + dt; t < endtime_lcl; t += dt){


        //snapshot code
        if(mod_cnt == snap_to_snap or mod_cnt == 0){
    
            mod_cnt = 0;
            progress = progress + 1;
            cout << 1.0*progress/(samplesize*snapshots_number) << endl;

            for(int x = 0; x < size; x++){
                psi_snapshots[x][snap_cnt] = psi[x];
            }

            for(int i = 0; i < size; i++){
                trace[snap_cnt] += conj(psi_snapshots[i][snap_cnt])*psi_snapshots[i][snap_cnt];
            }
            snap_cnt++;
        }
        mod_cnt++;


        EvolveStep(psinew, psi, p);
        
        tnorm[i] = 0.0;
        for (int j = 0; j < size; j++){
            psi[j] = psinew[j]; 
            tnorm[i] = tnorm[i] + pow(abs(psi[j]), 2); 
        }

        i++;



    } 
    //end of evolution

     
    //produce qubit density matrix
    complex <double> rho_qubit[transmon_levels][transmon_levels][snapshots_number];
    for(int i = 0; i < transmon_levels; i++){
        for(int j = 0; j < transmon_levels; j++){
            for(int k = 0; k < snapshots_number; k++){
                rho_qubit[i][j][k] = 0.0;
            }
        }
    }

    for(int i = 0; i < transmon_levels; i++){
        for(int j = 0; j < transmon_levels; j++){
            for(int k = 0; k < snapshots_number; k++){
                for(int n = 0; n < cavity_levels; n++){
                    rho_qubit[i][j][k] += conj(psi_snapshots[transmon_levels*n + i][k])*psi_snapshots[transmon_levels*n + j][k];
                }
            }
        }
    }

    //calculate Shannon entropy;
    /*
    complex <double> Trace[snapshots_number];
    complex <double> Det[snapshots_number];  
    complex <double> lambda1[snapshots_number];
    complex <double> lambda2[snapshots_number];
    for(int k=0; k < snapshots_number; k++){
        Trace[k]= rho_qubit[1][1][k] + rho_qubit[0][0][k];
        Det[k]= rho_qubit[1][1][k]*rho_qubit[0][0][k] - rho_qubit[1][0][k]*rho_qubit[0][1][k];
        lambda1[k]=(+Trace[k] + sqrt(Trace[k]*Trace[k]-4.*Det[k]))/2.;
        lambda2[k]=(+Trace[k] - sqrt(Trace[k]*Trace[k]-4.*Det[k]))/2.;
        S[k]=lambda1[k]*log(lambda1[k]) + lambda2[k]*log(lambda2[k]);
    }
    */
    /*
    //find sigmas for each snapshot
    for(int k=0; k<snapshots_number; k++){
        sigmaz[k]+=rho_qubit[1][1][k] - rho_qubit[0][0][k];
        sigmaplus[k]+=rho_qubit[1][0][k];
        sigmaminus[k]+=rho_qubit[0][1][k];
    }
    */

    //add to photon histogram
    for(int i = 0; i < (cavity_levels); i++){
        for(int k = 0; k < snapshots_number; k++){
            for( int j = 0; j < transmon_levels; j++){
                photon_histogram[i][k] += real(conj(psi_snapshots[transmon_levels*i + j][k])*psi_snapshots[transmon_levels*i + j][k]);
            }
        }
    }

    //add to transmon histogram
    for(int i = 0; i < transmon_levels; i++){
        for(int k = 0; k < snapshots_number; k++){
            transmon_histogram[i][k] += real(rho_qubit[i][i][k]);
        }
    }

    //find <n>
    for(int k=0; k<snapshots_number; k++) {
        for(int i=0; i < (cavity_levels); i++){
            for(int j = 0; j < transmon_levels; j++){
                photons[k] += i*real(conj(psi_snapshots[transmon_levels*i + j][k])*psi_snapshots[transmon_levels*i + j][k]);
            }
        }
    }

    //for ( int a = 0; a < glb_size; a++ ) {
    //    cout << setprecision(10) << psi_snapshots[a][0] << endl;
    //}

    //find <a> //fix this
    for(int k=0; k<snapshots_number; k++) {
        for(int i = 0; i < (cavity_levels) - 1; i++) {
            for(int j = 0; j < transmon_levels; j++) {
                averagealpha[k] += sqrt(i + 1)*(conj(psi_snapshots[transmon_levels*i + j][k])*psi_snapshots[transmon_levels*i + 2 + j][k]);
            }
        }
    }

    //find <a\sigma^{+}>, invalid for transmon
    /*
    for(int k = 0; k < snapshots_number; k++) {
        for(int i = 0; i < (cavity_levels) - 1; i++) {
            averagealphasigmaplus[k] += sqrt(i + 1)*conj(psi_snapshots[2*i + 1][k])*psi_snapshots[2*i + 2][k];
        }
    }
    */

    sample_cnt += 1;
	 
    #pragma omp critical
    {

    string path1 = "results/" + str + "/photons9.dat";
    ofstream myfile1 (path1.c_str(),ios_base::app | ios_base::out);
    myfile1 << scientific;
    myfile1 << setprecision(precision);
	for (int i=0; i<snapshots_number; i++) {
        myfile1 << real(photons[i]) << ',';
    }
    myfile1 << endl;
    myfile1.close();

/*  string path3 = "results/" + str + "/report.dat";    
    ofstream myfile3 (path3.c_str(),ios_base::app | ios_base::out);
    myfile3 << scientific;
    myfile3 << setprecision(precision);
    for (int i=0; i<bins; i++) {
        myfile3 << i*p.dt << '\t';
        myfile3 << tnorm[i] << '\t'; 
        myfile3 << endl;
    }
    myfile3.close(); 
    
   
    string path4 = "results/" + str + "/trace.dat";    
    ofstream myfile4 (path4.c_str(),ios_base::app | ios_base::out);
    myfile4 << scientific;
    myfile4 << setprecision(precision);
    for(int i=0; i<snapshots_number;i++){
        myfile4 << real(trace[i])/samplesize << '\t';
        myfile4 << endl;}  
    myfile4.close();  
    
    string path5 = "results/" + str + "/m_elements.dat";
    ofstream myfile5 (path5.c_str(),ios_base::app | ios_base::out);
    myfile5 << scientific;
    myfile5 << setprecision(precision);
    for(int i=0; i<snapshots_number;i++){
        myfile5 << m_element[i] << '\t';
        myfile5 << endl;
    }
    myfile5.close();*/ 

    string path6 = "results/" + str + "/qubitsigmaz9.dat"; 
    ofstream myfile6 (path6.c_str(),ios_base::app | ios_base::out);
    myfile6 << scientific;
    myfile6 << setprecision(precision);
    for( int i=0; i<snapshots_number; i++ ) {
        myfile6 << real(sigmaz[i]) << ',';
    }
    myfile6 << endl;
    
    myfile6.close();  

    string path7 = "results/" + str + "/qubitsigmaplusR9.dat";
    ofstream myfile7 (path7.c_str(),ios_base::app | ios_base::out);
    myfile7 << scientific;
    myfile7 << setprecision(precision);
    for( int i=0; i<snapshots_number; i++ ) {
        myfile7 << real(sigmaplus[i]) << ',';
    }
    myfile7 << endl;
 
    myfile7.close(); 

    string path8 = "results/" + str + "/qubitsigmaminusR9.dat"; 
    ofstream myfile8 (path8.c_str(),ios_base::app | ios_base::out);
    myfile8 << scientific;
    myfile8 << setprecision(precision);
    for( int i=0; i<snapshots_number; i++ ) {
        myfile8 << real(sigmaminus[i]) << ',';
    }
    myfile8 << endl;

    myfile8.close(); 

    string path9 = "results/" + str + "/qubitsigmaplusI9.dat";     
    ofstream myfile9 (path9.c_str(),ios_base::app | ios_base::out);
    myfile9 << scientific;
    myfile9 << setprecision(precision);
    for( int i=0; i<snapshots_number; i++ ) {
        myfile9 << imag(sigmaplus[i]) << ',';
    }
    myfile9 << endl;

    myfile9.close(); 

    string path10 = "results/" + str + "/qubitsigmaminusI9.dat";     
    ofstream myfile10 (path10.c_str(),ios_base::app | ios_base::out);
    myfile10 << scientific;
    myfile10 << setprecision(precision);
    for( int i=0; i<snapshots_number; i++ ) {
        myfile10 << imag(sigmaminus[i]) << ',';
    }
    myfile10 << endl;
    myfile10.close(); 
    
    string path11 = "results/" + str + "/ReaverageAlpha9.dat";     
    ofstream myfile11 (path11.c_str(),ios_base::app | ios_base::out);
    myfile11 << scientific;
    myfile11 << setprecision(precision);
    for(int i=0; i<snapshots_number;i++){
        myfile11 << real(averagealpha[i]) << ',';
    }
	myfile11 << endl;
    myfile11.close();     

    string path12 = "results/" + str + "/ImaverageAlpha9.dat";        
    ofstream myfile12 (path12.c_str(),ios_base::app | ios_base::out);
    myfile12 << scientific;
    myfile12 << setprecision(precision);
    for(int i=0; i<snapshots_number;i++){
        myfile12 << imag(averagealpha[i]) << ',';
    }
    myfile12 << endl;
    myfile12.close();

/*
    string path15 = "results/" + str + "/Ralphasigmaplus9.dat";  
    ofstream myfile15 (path15.c_str(),ios_base::app | ios_base::out);
    myfile15 << scientific;
    myfile15 << setprecision(precision);
    for(int i=0; i<snapshots_number;i++){
        myfile15 << real(averagealphasigmaplus[i]) << ',';
    }
    myfile15 << endl;
    myfile15.close(); 

    string path16 = "results/" + str + "/Ialphasigmaplus9.dat";  
    ofstream myfile16 (path16.c_str(),ios_base::app | ios_base::out);
    myfile16 << scientific;
    myfile16 << setprecision(precision);
    for(int i=0; i<snapshots_number;i++){
        myfile16 << imag(averagealphasigmaplus[i]) << ',';
    }
    myfile16 << endl;
    myfile16.close(); 

    string path17 = "results/" + str + "/RhocavityRe9.dat";  
    ofstream myfile17 (path17.c_str(),ios_base::app | ios_base::out);
    myfile17 << scientific;
    myfile17 << setprecision(precision);
    for(int k=0; k<snapshots_number;k++){
        for(int i=0;i<(p.size/2);i++){
            for(int j=0;j<(p.size/2);j++){
                myfile17 << real(rho_cavity[i][j][k]) << '\t';
            }
        myfile17 <<";"<< endl;
        }
    myfile17 <<";"<< endl;
    }
    myfile17.close(); 

    string path18 = "results/" + str + "/RhocavityIm9.dat";  
    ofstream myfile18 (path18.c_str(),ios_base::app | ios_base::out);
    myfile18 << scientific;
    myfile18 << setprecision(precision);
    for(int k=0; k<snapshots_number;k++){
        for(int i=0;i<(p.size/2);i++){
            for(int j=0;j<(p.size/2);j++){
                myfile18 << imag(rho_cavity[i][j][k]) << '\t';
            }
        myfile18 <<";"<< endl;
        }
    myfile18 <<";"<< endl;
    }
    myfile18.close();*/

    //not used for transmons
    /*
    string path19 = "results/" + str + "/Entropy.dat";  
    ofstream myfile19 (path19.c_str(),ios_base::app | ios_base::out);
    myfile19 << scientific;
    myfile19 << setprecision(precision);
    for(int i=0; i<snapshots_number; i++){
        myfile19 << real(S[i]) << ',';
    }
	myfile19 << endl;
    myfile19.close();
    */

    }
    
}



int main(){

    using namespace std::chrono;
    system_clock::time_point tp = system_clock::now();
    system_clock::duration dtn = tp.time_since_epoch();

    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,80,"%Y.%m.%d--%H-%M-%S",timeinfo);
    string str(buffer);
    str = folder + '/' + str;
    string subfolder = "mkdir results/" + folder;
    system(subfolder.c_str());
    string results_folder = "mkdir results/" + str;
    system(results_folder.c_str());

    for(int i = 0; i < cavity_levels; i++){
        for(int k = 0; k < snapshots_number; k++){
            photon_histogram[i][k] = 0.0;
        }
    }

    for(int i = 0; i < transmon_levels; i++){
        for(int k = 0; k < snapshots_number; k++){
            transmon_histogram[i][k] = 0.0;
        }
    }

    ofstream parameters_file;
    parameters_file << scientific;
    parameters_file << setprecision(precision);
    ofstream settings_file;
    settings_file << scientific;
    settings_file << setprecision(precision);

    myfunc_params p;
    p.kappa = kappa;
    p.gamma = gamma_c;
    p.g = g;
    p.dt = dt;
    p.size = glb_size;
    p.drive = epsilon;
    p.omega_c = omega_c;
    p.omega_q = omega_q;
    p.omega_d = omega_d;
    p.dc = dc;
    p.dq = dq;
    p.chi = chi;

    const gsl_rng_type * Trng;

    gsl_rng_env_setup();
    Trng = gsl_rng_default;
    r = gsl_rng_alloc (Trng);
    gsl_rng_set(r, dtn.count());
    //gsl_rng_set(r, 19);
    char filename[100];

    string parameters_path = "results/" + str + "/parameters.txt";
    parameters_file.open(parameters_path.c_str());
    parameters_file << "Drive_frequency " << p.omega_d << endl;
    parameters_file << "Qubit_frequency " << p.omega_q << endl;
    parameters_file << "Cavity_frequency " << p.omega_c << endl;
    parameters_file << "Coupling " << p.g << endl;
    parameters_file << "Kappa " << p.kappa << endl;
    parameters_file << "Gamma " << p.gamma << endl;
    parameters_file << "Endtime " << endtime << endl;
    parameters_file << "Dt " << p.dt << endl;
    parameters_file << "Bins " << bins << endl;
    parameters_file << "Size " << p.size << endl;
    parameters_file << "Snap_to_snap " << snap_to_snap << endl;
    parameters_file << "Snapshots_number " << snapshots_number << endl;
    parameters_file << "Epsilon " << epsilon << endl;
    parameters_file << "Theta " << theta << endl;
    parameters_file << "Phi " << phi << endl;
    parameters_file.close();

    string settings_path = "results/" + str + "/settings.cfg";
    settings_file.open(settings_path.c_str());
    settings_file << "folder=" << folder << endl;
    settings_file << "endtime=" << endtime << endl;
    settings_file << "transmon_levels=" << transmon_levels << endl;
    settings_file << "cavity_levels=" << cavity_levels << endl;
    settings_file << "snap_to_snap=" << snap_to_snap << endl;
    settings_file << "bins=" << bins << endl;
    settings_file << "snapshots_number=" << snapshots_number << endl;
    settings_file << "kappa=" << kappa << endl;
    settings_file << "gamma_c=" << gamma_c << endl;
    settings_file << "g=" << g << endl;
    settings_file << "epsilon=" << epsilon << endl;
    settings_file << "omega_c=" << omega_c << endl;
    settings_file << "omega_q=" << omega_q << endl;
    settings_file << "omega_d=" << omega_d << endl;
    settings_file << "theta=" << theta << endl;
    settings_file << "phi=" << phi << endl;
    settings_file << "samplesize=" << samplesize << endl;
    settings_file << "photon_number=" << photon_number << endl;
    settings_file << "chi=" << chi << endl;
    settings_file.close();

    #pragma omp parallel for num_threads(50)
    for (int i = 0; i < samplesize; i++) //samplesize is the number of trajectories
    {
        TrajectoryAdd(&p,endtime,i,str);
    }

    for(int i = 0; i < cavity_levels; i++){
        for(int k = 0; k < snapshots_number; k++){
            photon_histogram[i][k] = photon_histogram[i][k]/samplesize;
        }
    }

    for(int i = 0; i < transmon_levels; i++){
        for(int k = 0; k < snapshots_number; k++){
            transmon_histogram[i][k] = transmon_histogram[i][k]/samplesize;
        }
    }

    string path20 = "results/" + str + "/photon_histogram.dat";  
    ofstream myfile20 (path20.c_str(),ios_base::app | ios_base::out);
    myfile20 << scientific;
    myfile20 << setprecision(precision);
    for(int i = 0; i < cavity_levels; i++){
        for(int k = 0; k < snapshots_number; k++){
            myfile20 << photon_histogram[i][k] << ',';
        }
        myfile20 << endl;
    }
    myfile20.close();

    string path21 = "results/" + str + "/transmon_occupations.dat";
    ofstream myfile21 (path21.c_str(),ios_base::app | ios_base::out);
    myfile21 << scientific;
    myfile21 << setprecision(precision);
    for(int i = 0; i < transmon_levels; i++){
        for(int k = 0; k < snapshots_number; k++){
            myfile21 << transmon_histogram[i][k] << ',';
        }
        myfile21 << endl;
    }
    myfile21.close(); 


   return EXIT_SUCCESS;
}
