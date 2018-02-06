#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <cstring>
#include <malloc.h>
#include <algorithm>
#include <omp.h>

using namespace std;
#define PI (atan(1)*4)

//-----------------------------------------------------------------------------
//		Parameters with fixed value 
//-----------------------------------------------------------------------------
#define G_Na 120
#define E_Na 50
#define G_K 36
#define E_K -77
#define G_L 0.3
#define E_L -54.387	
#define C 1
#define V_G_E 0
#define V_G_I -80
#define Sigma_r_E 0.5
#define Sigma_r_I 0.5
#define Sigma_d_E 3.0
#define Sigma_d_I 7.0

double Sigma_ratio_E = Sigma_d_E*Sigma_r_E / (Sigma_d_E - Sigma_r_E);
double Sigma_ratio_I = Sigma_d_I*Sigma_r_I / (Sigma_d_I - Sigma_r_I);

#define V_th -50.0
#define Epsilon 1e-10
#define T_ref 3.5

//-----------------------------------------------------------------------------
//		Variables for input parameters
//-----------------------------------------------------------------------------
int N, NE, NI;                // Number of neurons(total,excitatory,inhibirory)
double T_Max, T_Step_Large, T_Step_Small;         // Total time & time step
double S[4];                  // Coupling strength (E->E,E->I,I->E,I->I)
int I_CONST;                  // electrode current constant  
double I_const_input;         // constant input current


double Nu,f;               // Feedforward Poisson rate and strength(E,I)  
double P_c;                   // Connect probability
int random_S, random_Nu;					
int IntOrd=1;					  //interpolation_order used in library method,1--linear,2--quadratic,3--cubic	
int method;					  //0-RK,1-Ad,2-Lib,3-ETD_RK,4-ETD	
int Regular_method, Lib_method, Adaptive_method,ETDRK_method;  
int Lyapunov;                     // compute largest lyapunov exponnet
int Power_spectrum = 0;				  // record v for power spectrum	
int Estimate_RK4_call = 0;					// compare the efficiency
int record_data[2];				// save data or not
char input_parameters[200];     //file for the parameters
char file[200],file1[200];				  // Record data path & Library path
int RecordFP = 0;					// record fire pattern 0101000101

//-----------------------------------------------------------------------------
//		Netwrok Information
//-----------------------------------------------------------------------------
double **Connect_Matrix;         // Connect matrix
double **CS;					 // Coupling strength matrix
struct neuron 
{
	double t,Nu;
	double v,dv,m,n,h,G_se,G_sse,G_si,G_ssi,G_f,G_ff;
	double I_input;
	double last_fire_time;
	int fire_num;
	int if_fired;
	double *Poisson_input_time;
	int Poisson_input_num;
	long seed;
	double wait_strength_E, wait_strength_I;
	double last_hit[4]; // I,m,h,n last fired with library method

	double F[4][9]; //ETD method, history, t,v,m,h,n,dv,dm,dh,dn
	int id_F;

	int state;     //1--neu,0--neu_old
};
struct neuron *neu, *neu_old;


#define MIN(a,b)  ((a)<(b)?(a):(b))
#define MAX(a,b)  ((a)>(b)?(a):(b))
//-----------------------------------------------------------------------------
//		Record firing time and voltage
//-----------------------------------------------------------------------------
FILE *FP,*FP1, *FP_FFTW, *FP_fire_pattern;
FILE *ffp; // for test
FILE *fp_v;  // usd when build the library, trace of V 
int library_v_trace = 0;
//-----------------------------------------------------------------------------
//		Library. Lib_input include I_input,m,h,n,V,m,h,n
//      input I_input,m,n,h, we can get V,m,h,n after evolving Tref
//      Lib_length, resolution record the data length and unique resolution
//-----------------------------------------------------------------------------
double **Lib_data,**Lib_unique;  
int Lib_length, *Lib_resolution;
double **Lib_v;

// the number of calling RK4 calculation & synchronize num (during T_ref)
double Call_num = 0, syn_num = 0, count_num = 0;
