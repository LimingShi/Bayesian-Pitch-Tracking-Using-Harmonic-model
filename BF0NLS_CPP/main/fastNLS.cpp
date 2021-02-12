#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include  <stdio.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "../src/tools.hpp"
#include "gtest/gtest.h"
//#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

#include "../src/single_pitch.hpp"

#ifdef DOUBLE
#define EPS 1e-10 // pure floating
#define EPS2 1e-8 //iterative algs
#else
#define EPS 1e-1
#define EPS2 1e-1
#endif



double* readInput(char *filename, int *numSample, int *fs)
{
    hid_t file;
    herr_t status;
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    double *nData=vector(1);
    H5LTread_dataset(file, "/DL1", H5T_NATIVE_DOUBLE, nData);
    *numSample=nData[0];
    
    double * x = vector(*numSample);
    H5LTread_dataset(file, "/DS1", H5T_NATIVE_DOUBLE, x);
    
    double *sr=vector(1);
    H5LTread_dataset(file, "/SR1", H5T_NATIVE_DOUBLE, sr);
    *fs=sr[0];
    
    status = H5Fclose(file);
    
    double sumE=0;
    for( int ii = 1 ; ii <= nData[0] ; ++ii){
        sumE += (x[ii-1])*(x[ii-1]);
    }
    sumE = sqrt(sumE/nData[0]);
    double scale = sqrt(3.1623e-5)/sumE; // scale to -45-dB loudness level
    //     double scale=1;
    
    for (int ii = 1 ; ii <= nData[0] ; ++ii)
        x[ii-1] = scale*x[ii-1];
    return x;
    
}


void process(char *filename, char *out_put_filename)
{
    int sigLength; int F;
    
    double *x=readInput(filename, &sigLength, &F);
    
    
    
    
    
    
    int N = (int) round(25*F/1000/2)*2;
    
    int L = (int) 10;
    
    int shift_vaue= (int) round(10*F/1000);
    //    zero padding
    x=x-N/2;sigLength=sigLength+N/2;
    int nData = sigLength;
    for (int ii =0; ii<N/2;++ii)
    {
        x[ii]=0.0;
    }
    
    
    
    
    int nframes= (int) floor(((double)nData-N)/shift_vaue)+1;
    
    
    int Mp = F/2+1;
    char label[50];
    
    
    
    double pitch_bounds[2];
    pitch_bounds[0]=70.0/F; pitch_bounds[1]=400.0/F;
    //   H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
#ifdef DOUBLE
    double * pb = pitch_bounds;
#else
    float pb[2];
    pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
#endif
    
    double trans_std_var = 2.0/F;
    //     double trans_std_var = 1e9/(double)F;
    
    FTYPE voicing_prob=0.7;
    
    /* construct object and that will compute gamma1 and gamma2  */
    single_pitch * sp = new single_pitch(L, N, pb,trans_std_var, voicing_prob);
    
    
    double * xp= vector(N);
    FTYPE * omega_0h;
    int estModelOrder;
    
    std::ofstream myfile1 (out_put_filename);
    
     myfile1 << "time "<<"F0 "<<"order "<<"Voicing "<<"F0_1 "<<"order_1"<<"\n";
    
    for( int ii = 1 ; ii <= nframes ; ++ii){
        
        CvmdCopy(N, x, xp);
        
        omega_0h = sp->est_fast(xp);
        

        myfile1 << (ii-1)*10/1000.0;
        myfile1 << " ";
        
        myfile1 << omega_0h[0]*F/2.0/M_PI;
        myfile1 << " ";
        myfile1 << omega_0h[1];
        myfile1 << " ";
        myfile1 << omega_0h[2];
        
        myfile1 << " ";
        if (omega_0h[2]>0.5)
        {
            myfile1 << omega_0h[0]*F/2.0/M_PI;
        }
        else
        {
            myfile1 << 0.0/0.0;
        }
        myfile1 << " ";
        
        if (omega_0h[2]>0.5)
        {
            myfile1 << omega_0h[1];
        }
        else
        {
            myfile1 << 0.0/0.0;
        }
        
        myfile1 << "\n";
        
        x=x+shift_vaue;
    }
    myfile1.close();
    delete sp;
}

int main(int argc, char **argv){
    if (argc!=3)
    {
        printf("Usage: fastNLS input output\n");
        return 1;
    }
    //     int F=0;
    //
    //
    //     char *s= argv[1];
    //
    //     sscanf(s, "%d", &F);
    
    //     int F=argv[1];;
    char *filename=argv[1];
    char *outfn=argv[2];
    
    
    process(filename,outfn);
    
    return 0;
    
}

