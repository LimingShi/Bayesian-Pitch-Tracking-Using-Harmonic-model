#ifndef __SINGLE_PITCH_H__
#define __SINGLE_PITCH_H__


#include <math.h>
#include <complex.h>
#include <limits>
#include <iostream>

#ifdef MKL
#include "fftw3.h"
#else
#include <fftw3.h>
#endif

#include "tools.hpp"
#include "th.hpp"
#include "vector.hpp"
//#include <Eigen/Dense>
//#include <Eigen/Core>

#ifdef REGULARIZE
#define T0_REG 1e-8
#else
#define T0_REG 0.0
#endif

#define DIM int

//using namespace Eigen;

class single_pitch{
public:
  single_pitch(int maxModelOrder, int nData, 
                 FTYPE * pitchBounds, double trans_std_var,FTYPE voicing_prob);

  ~single_pitch();

  void compute(FTYPE * x);
  void update_gamma(int iOrder, int nPitches, int nPitchesOld, 
                    int nPitchesOldOld, FTYPE * phi, 
                    FTYPE * psi, FTYPE * alpha, FTYPE * lambda, 
                    FTYPE * mu, FTYPE * Gamma);

  void update_gamma_p(int iOrder, int nPitches, int nPitchesOld, 
                      int nPitchesOldOld, FTYPE * alpha, FTYPE * Gamma);

  void update_ls_sol(int iOrder, int nPitches, int nPitchesOld, 
                     bool add, FTYPE * dftData, 
                     FTYPE * Gamma, FTYPE * lsSol);

  void compute_cc();
  
//  added by Liming Shi
  inline double squared(double x);
    double normpdf(double x, double u, double s);
    double log_sum_exp(double arr[], int count);
    double log_sumsum_exp(FTYPE **arr, int row_num,int col_num, double A_num);
//    double ramp(double x);
    double norm_prob(FTYPE **arr, int row_num,int col_num, double A_num,double scale_par);
    int norm_prob2(FTYPE **arr, int row_num, int col_num, double scale_par,FTYPE **arrtarget);
    double log_sum(double arr[], int count);
    double sum_vec(double arr[], int count);
    //
  
  FTYPE compute_obj(FTYPE omega, FTYPE * x, int iOrder);
  FTYPE compute_obj(FTYPE omega, FTYPE * x, int iOrder, FTYPE * ac, FTYPE * as);

  FTYPE * refine(FTYPE * x);
  FTYPE * refine(FTYPE * x, FTYPE eps);
  FTYPE refine_single(FTYPE *x, int iOrder, FTYPE eps);
  void compute_max_on_grid(void);
    
  FTYPE golden(FTYPE * x, int iOrder, FTYPE omega_L, FTYPE omega_U, FTYPE eps);

//   int model_order_selection();
  int model_order_selection();

  FTYPE est(FTYPE * x);
  FTYPE est(FTYPE * x, FTYPE lnBFZeroOrder, FTYPE eps);
  FTYPE * est_fast(FTYPE * x);

  FTYPE max_obj(int iOrder);
  FTYPE argmax_obj(int iOrder);

  FTYPE * Gamma1;
  FTYPE * Gamma2;
  FTYPE ** costFunctionMatrix;
  FTYPE * omega_0h;
  FTYPE * lnBF;

  int nPitches(int iOrder){ return nPitchesAll[iOrder-1];};
  int modelOrder(){return estOrder;};

private:
  
  int maxFftIndexOld;
//int M, N, L, F;
  int maxModelOrder, nData, nFftGrid;
  int estOrder;

  int t;
  int k;
  int M;
  int Mp;
  int minFftIndex;

  FTYPE pitchBounds[2];

  int * nPitchesAll;

  FTYPE ** crossCorrelationVectors;
  FTYPE ** tf;
  FTYPE ** hf;

  FTYPE * dftData1;
  FTYPE * dftData2;
  FTYPE * x;
  FTYPE * y;
  FTYPE * yt;

  FTYPE * t1;
  FTYPE * t2;
  FTYPE * t3;
  FTYPE * t4;
  FTYPE * tp;
 
  FTYPE * fftShiftVector;
  FTYPE * lsSol1;
  FTYPE * lsSol2;

  FTYPE * range;

  FTYPE * xp;
  FTYPE * dftData;
  FFTW_PLAN p1;
//   added by Liming Shi
  FTYPE * logpi;
    FTYPE ** A;
    FTYPE ** B;
    FTYPE ** C;
    FTYPE * temp_normalization;
    FTYPE ** scaled_alpha_buffer;
    FTYPE ** scaled_alpha_buffer2;
    FTYPE unvoicing_scaled_alpha_buffer;
    FTYPE ** bar_alpha;
    FTYPE unvoicing_scaled_alpha;
//    FTYPE bar_alpha;
    FTYPE * state_prior;
    FTYPE * logModelPrior;
    FTYPE * returned_value;
};


#ifdef LIB
// The class is wrapped using these C-style functions
extern "C"{

  single_pitch * single_pitch_new(int maxModelOrder, int nFftGrid, 
                                  int nData, FTYPE * pitchBounds){
    return new single_pitch(maxModelOrder, nFftGrid, 
                            nData, pitchBounds);
  }
  
  FTYPE single_pitch_est(single_pitch * sp, FTYPE * x, FTYPE lnBFZeroOrder, FTYPE eps){
    return sp->est(x, lnBFZeroOrder, eps);
  }

  FTYPE single_pitch_est_fast(single_pitch * sp, FTYPE * x, FTYPE lnBFZeroOrder, FTYPE eps){
    return sp->est_fast(x, lnBFZeroOrder, eps);
  }
  
  void single_pitch_del(single_pitch * sp){delete sp;}

  int single_pitch_model_order(single_pitch * sp){return sp->modelOrder();}
}
#endif

#endif
