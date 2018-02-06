
/* Hodgkin-Huxley model: Network (Library and Regular)*/
//-----------------------------------------------------------------------------
//		Comments
//-----------------------------------------------------------------------------


#include "Def.h"
#include "Gate_variables.h"
#include "Random.h"
#include "Read_parameters.h"
#include "Initialization.h"
#include "Find_cubic_hermite_root.h"
#include "Runge_Kutta4.h"
#include "ETDRK.h"
#include "test_func.h"
#include "Run_model.h"
#include "Largest_Lyapunov.h"
#include "Total_call_RK4.h"
#include "Delete.h"





int main(int argc,char **argv)
{
	long seed, seed0, seed1, seed2;
	clock_t t0, t1;	 
	char str[200];
	double MLE;
	double mean_fire_rate;

	strcpy(input_parameters, argv[1]);   
	Read_parameters(seed, seed1);

	/////////////////////	
	out_put_filename();
	seed0 = seed;    // Create connect matrix
	seed2 = seed1;  // Initialization & Poisson
	Initialization(seed0, seed2);

	t0 = clock();
	if (Lyapunov)
		MLE = Largest_Lyapunov(seed2, 1, T_Step_Large);
	else
		Run_model();


	int total_fire_num[2] = { 0 };
	for (int i = 0; i < N; i++)
		total_fire_num[i < NE ? 0 : 1] += neu[i].fire_num;

	
	mean_fire_rate = (total_fire_num[0] + total_fire_num[1]) / T_Max * 1000 / N; //(Hz)
	printf("mean firing rate = %0.2f(Hz)\n", mean_fire_rate);

	if (Lyapunov)
	{
		printf("The largest Lyapunov exponent is %f\n", MLE);
		return MLE;
	}
	//t1 = clock();
	if (Estimate_RK4_call)
	{
		double s = new_total_call_RK4(mean_fire_rate / 1000);
		printf("num1=%0.0f num2=%0.0f %f\n", Call_num, s, s / Call_num);
	}

	//printf("Total time = %0.3fs \n\n", double(t1 - t0) / CLOCKS_PER_SEC);
	Delete();
	
}


