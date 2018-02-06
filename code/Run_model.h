/////////////// evolve for one single time step from t to t+dt
void evolve_model_with_correct_timestep(int n, struct neuron &a, double t, double dt)
{	
	if (dt < 0)
	{
		printf("Warning! dt in evolve_model_with_correct_timestep < 0 !! dt = %e\n", dt);
		printf("n=%d t=%f state=%d\n", n, t,a.state);
		exit(0);
	}

	if (Regular_method || Adaptive_method)
	{
		Update_RK4(n, a, t, dt);
	}
	else if (Lib_method)
	{
		double local_t = a.last_fire_time + T_ref;
		if (local_t <= t)
		{
			Update_RK4(n, a, t, dt);
		}
		else if (local_t > t + dt)
			Update_neu_G(n, a, t, dt);
		else  // t < local_t < t+dt case
		{
			Update_neu_G(n, a, t, local_t - t);
			Update_RK4(n, a, local_t, t + dt - local_t);
		}
	}
	else
	{
		Update_ETDRK4(n, a, t, dt);
	}
	
}


/////////////// evolve from t to t+dt 
///////////// there may be times nodes in [t,t+dt], from Poisson and Adaptive time nodes
void Check_update_conti_or_Poisson_input(int n, struct neuron &a, double t, double dt)
{
	if (Nu < Epsilon)
		evolve_model_with_correct_timestep(n, a, t, dt);
	else
	{
		double t1 = t;
		for (int i = 0; i < neu[n].Poisson_input_num; i++)
		{
			if (neu[n].Poisson_input_time[i] >= t1 && neu[n].Poisson_input_time[i] < t + dt)
			{
				evolve_model_with_correct_timestep(n, a, t1, neu[n].Poisson_input_time[i] - t1);
				t1 = neu[n].Poisson_input_time[i];
				a.G_ff += f;
			}
		}
		evolve_model_with_correct_timestep(n, a, t1, t + dt - t1);
	}
}


void evolve_with_multi_time_nodes(int n, struct neuron &a, double t, double dt)
{
	if (Regular_method || Lib_method || ETDRK_method || t >= a.last_fire_time + T_ref)  // no time nodes from T_Step_Small
		Check_update_conti_or_Poisson_input(n, a, t, dt);
	else
	{
		double DT, local_t;
		
		local_t = a.last_fire_time + T_ref;

		if (dt > local_t - t)
		{
			DT = local_t - t;
			int step = int(DT / T_Step_Small);

			for (int k = 0; k < step; k++)
				Check_update_conti_or_Poisson_input(n, a, t + k*T_Step_Small, T_Step_Small);
			Check_update_conti_or_Poisson_input(n, a, t + step*T_Step_Small, local_t - t - step*T_Step_Small);   ///To local_t
			Check_update_conti_or_Poisson_input(n, a, local_t, t + dt - local_t);		
		}
		else
		{
			DT = dt;
			int step = int(DT / T_Step_Small);

			for (int k = 0; k < step; k++)
				Check_update_conti_or_Poisson_input(n, a, t + k*T_Step_Small, T_Step_Small);
			Check_update_conti_or_Poisson_input(n, a, t + step*T_Step_Small, dt - step*T_Step_Small);
		}

	}
}

void Exchange(struct neuron &a, struct neuron b)  // a <-- b
{
	a.t = b.t;
	a.Nu = b.Nu;
	a.v = b.v;
	a.dv = b.dv;
	a.m = b.m;
	a.h = b.h;
	a.n = b.n;
	a.G_se = b.G_se;
	a.G_sse = b.G_sse;
	a.G_si = b.G_si;
	a.G_ssi = b.G_ssi;
	a.I_input = b.I_input;

	a.last_fire_time = b.last_fire_time;
	a.G_f = b.G_f;
	a.G_ff = b.G_ff;

	a.if_fired = 0;
}


//-----------------------------------------------------------------------------
//		Evolve single neuron from last hit time t to t+T_ref with library method (v,m,h,n)
//      Multi interpolation with interpolation_order (1--linear,2--quadratic)
//		Multi-interpolation, input dimension, outout dimension
//-----------------------------------------------------------------------------

void Update_Lib_way(struct neuron &neu_a, int interpolation_order, int input_dim, int output_dim)
{
	// I_input,m,h,n ____ 4-digit number like decimalism but with each lenth l[4] not 10    
	// l--length in each I,m,h,n. this is actually Lib_resolution
	// a--each minimin id for interpolation in 4-D. 
	// b--each id for interpolation in 4-D. scan id
	// c--each change in id id for interpolation in 4-D. scan change

	int *l;  
	double *input, *output; //I,m,h,n & v,m,h,n
	int *a,*b, *c;

	input = new double[input_dim]; output = new double[output_dim];
	l = new int[input_dim]; a = new int[input_dim];
	b = new int[input_dim]; c = new int[input_dim];

	for (int i = 0; i < input_dim; i++)
		output[i] = 0;

	for (int i = 0; i < input_dim; i++)
		l[i] = Lib_resolution[i];

	input[0] = neu_a.I_input;
	input[1] = neu_a.m;
	input[2] = neu_a.h;
	input[3] = neu_a.n;

	for (int i = 0; i < input_dim; i++)
	{
		a[i] = int((input[i] - Lib_unique[i][0]) / (Lib_unique[i][1] - Lib_unique[i][0]) + 0.1);
		if (a[i] < 0)
		{
			printf("Warning! a%d=%d\n", i, a[i]);
			a[i] = 0;
		}
		else if (a[i] > l[i] - interpolation_order - 1)
		{
			printf("Warning! a%d=%d max_length=%d t=%0.3f\n ", i, a[i], Lib_resolution[i], neu_a.t);
			a[i] = l[i] - interpolation_order - 1;
		}
	}


	for (int i = 0; i < int(pow((interpolation_order + 1), input_dim) + 0.05); i++)
	{
		int k;
		double s = 1;
		k = i;
		for (int j = 0; j < input_dim; j++)
		{
			c[j] = k % (interpolation_order + 1);
			k /= (interpolation_order + 1);
			for (int ll = 0; ll < (interpolation_order + 1); ll++)
			{
				if (ll == c[j])
					continue;
				else
					s *= (input[input_dim - 1 - j] - Lib_unique[input_dim - 1 - j][a[input_dim - 1 - j] + ll]) /
					(Lib_unique[input_dim - 1 - j][a[input_dim - 1 - j] + c[j]] - Lib_unique[input_dim - 1 - j][a[input_dim - 1 - j] + ll]);
			}
		}
		for (int j = 0; j < input_dim; j++)
			b[j] = a[j] + c[input_dim - 1 - j];

		int m = b[0];
		for (int j = 1; j < input_dim; j++)
			m = m*l[j] + b[j];

		for (int j = 0; j < output_dim; j++)
			output[j] += Lib_data[j + input_dim][m] * s;
	}


	neu_a.v = output[0];
	neu_a.m = output[1];
	neu_a.h = output[2];
	neu_a.n = output[3];
	delete[]l, delete[]a, delete[]input, delete[]output, delete[]b, delete[]c;

}

//-----------------------------------------------------------------------------
//		To record v trace in Power spectrum in Lib_method
//-----------------------------------------------------------------------------

void Find_v_in_lib_method(struct neuron neu_a, int interpolation_order, int input_dim, double &output, double dt)
{
	// I_input,m,h,n ____ 4-digit number like decimalism but with each lenth l[4] not 10    
	// l--length in each I,m,h,n. this is actually Lib_resolution
	// a--each minimin id for interpolation in 4-D. 
	// b--each id for interpolation in 4-D. scan id
	// c--each change in id id for interpolation in 4-D. scan change

	int t_id = int(dt * 256);  // time_step = 2^-8, 0 <= dt <= T_ref
	double lambda_t = dt * 256 - t_id;

	int *l;
	double *input; //I,m,h,n
	int *a, *b, *c;


	output = 0;

	input = new double[input_dim];
	l = new int[input_dim]; a = new int[input_dim];
	b = new int[input_dim]; c = new int[input_dim];


	for (int i = 0; i < input_dim; i++)
	{
		l[i] = Lib_resolution[i];
		input[i] = neu_a.last_hit[i];
	}

	for (int i = 0; i < input_dim; i++)
	{
		a[i] = int((input[i] - Lib_unique[i][0]) / (Lib_unique[i][1] - Lib_unique[i][0]) + 0.1);
		if (a[i] < 0)
		{
			printf("Warning! a%d=%d\n", i, a[i]);
			a[i] = 0;
		}
		else if (a[i] > l[i] - interpolation_order - 1)
		{
			printf("Warning! a%d=%d max_length=%d t=%0.3f\n", i, a[i], Lib_resolution[i], neu_a.t);
			a[i] = l[i] - interpolation_order - 1;
		}
	}


	for (int i = 0; i < int(pow((interpolation_order + 1), input_dim) + 0.05); i++)
	{
		int k;
		double s = 1;
		k = i;
		for (int j = 0; j < input_dim; j++)
		{
			c[j] = k % (interpolation_order + 1);
			k /= (interpolation_order + 1);
			for (int ll = 0; ll < (interpolation_order + 1); ll++)
			{
				if (ll == c[j])
					continue;
				else
					s *= (input[input_dim - 1 - j] - Lib_unique[input_dim - 1 - j][a[input_dim - 1 - j] + ll]) /
					(Lib_unique[input_dim - 1 - j][a[input_dim - 1 - j] + c[j]] - Lib_unique[input_dim - 1 - j][a[input_dim - 1 - j] + ll]);
			}
		}

		for (int j = 0; j < input_dim; j++)
			b[j] = a[j] + c[input_dim - 1 - j];

		int m = b[0];
		for (int j = 1; j < input_dim; j++)
			m = m*l[j] + b[j];

		double find_v;
		find_v = Lib_v[t_id][m] * (1 - lambda_t) + Lib_v[t_id + 1][m] * lambda_t;
		output += find_v * s;

	}
	delete[]a, delete[]l, delete[]input, delete[]b, delete[]c;
}


//-----------------------------------------------------------------------------
//		find next firing time in [t, t1] for single neuron,  if found, recorded by a_old 
//      (to find the first spike time. i.e. spike-spike corrections)
//-----------------------------------------------------------------------------

void find_next_fire_time_for_single_neu(int id, struct neuron a, struct neuron &a_old, double t, double t1)
{
	Exchange(a_old, a);
	a_old.if_fired = 0;

	evolve_with_multi_time_nodes(id, a_old, t, t1 - t);

}


//-----------------------------------------------------------------------------
//		Update all neurons from neu.t --> t and chech if fires during [t, t1]
//-----------------------------------------------------------------------------

void update_all_neu(struct neuron *a, struct neuron *a_old, double t, double t1)
{
	#pragma omp parallel for    
	for (int id = 0; id < N; id++)
	{
		a[id].if_fired = 0;
		if ((a[id].wait_strength_E == 0) && (a[id].wait_strength_I == 0))  
			continue;

		evolve_with_multi_time_nodes(id, a[id], a[id].t, t - a[id].t);

		a[id].G_sse += a[id].wait_strength_E;
		a[id].G_ssi += a[id].wait_strength_I;

		a[id].wait_strength_E = 0;
		a[id].wait_strength_I = 0;	

		if (!a[id].if_fired)
			find_next_fire_time_for_single_neu(id, a[id], a_old[id], t, t1);
	}

	int continue_update = 0;
	for (int id = 0; id < N; id++)
	{
		if (a[id].if_fired)
		{
			continue_update = 1;
			if (record_data[0])
			{
				fwrite(&t, sizeof(double), 1, FP);
				double s = id;
				fwrite(&s, sizeof(double), 1, FP);
			}
			a[id].fire_num++;
			a[id].last_fire_time = t;
			if (a[id].v < V_th)
				a[id].v = V_th;

			for (int i = 0; i < N; i++)
			{
				if (i == id)
					continue;
				if (id < NE && i < NE)
					a[i].wait_strength_E += CS[id][i];
				else if (id < NE && i >= NE)
					a[i].wait_strength_E += CS[id][i];
				else if (id >= NE && i < NE)
					a[i].wait_strength_I += CS[id][i];
				else
					a[i].wait_strength_I += CS[id][i];
			}

			if (Power_spectrum && Lib_method)
			{
				a[id].last_hit[0] = a[id].I_input;
				a[id].last_hit[1] = a[id].m;
				a[id].last_hit[2] = a[id].h;
				a[id].last_hit[3] = a[id].n;
			}

			if (Estimate_RK4_call && a[id].t > 10)
			{
				for (int i = 0; i < N; i++)
				{
					if (Connect_Matrix[id][i] == 1)
					{
						count_num++;
						if (a[id].t - a[i].last_fire_time <= T_ref && id != i)
							syn_num++;
					}
				}
			}


			if (Lib_method)
				Update_Lib_way(a[id], IntOrd, 4, 4);     
			find_next_fire_time_for_single_neu(id, a[id], a_old[id], t, t1);
			
			if(Lib_method && a_old[id].if_fired)
				a_old[id].last_fire_time += T_ref;
		}
	}
	if(continue_update)
		update_all_neu(a, a_old, t, t1);
}

//-----------------------------------------------------------------------------
//		Event driven
//-----------------------------------------------------------------------------

void evolve_model_with_initial_timestep(struct neuron *a, struct neuron *a_old, double t, double dt)
{
	#pragma omp parallel for    
	for (int i = 0; i < N; i++)        // parallel
		find_next_fire_time_for_single_neu(i, a[i], a_old[i], t, t + dt);

	while (1)
	{
		double first_fire_time = t + dt;
		int first_fire_neu = -1;

		for (int i = 0; i < N; i++)
		{
			if (a_old[i].if_fired && a_old[i].last_fire_time < first_fire_time)
			{
				first_fire_time = a_old[i].last_fire_time;
				first_fire_neu = i; 
			}
		}

		if (first_fire_neu > -1)
		{
			int id = first_fire_neu;

			evolve_with_multi_time_nodes(id, a[id], a[id].t, first_fire_time - a[id].t);

			if (record_data[0])
			{
				fwrite(&first_fire_time, sizeof(double), 1, FP);
				double s = id;
				fwrite(&s, sizeof(double), 1, FP);
			}
			a[id].fire_num++;
			a[id].last_fire_time = first_fire_time;
			if (a[id].v < V_th)
				a[id].v = V_th;

			for (int i = 0; i < N; i++)
			{
				if (i == id)
					continue;
				if (id < NE && i < NE)
					a[i].wait_strength_E += CS[id][i];
				else if (id < NE && i >= NE)
					a[i].wait_strength_E += CS[id][i];
				else if (id >= NE && i < NE)
					a[i].wait_strength_I += CS[id][i];
				else
					a[i].wait_strength_I += CS[id][i];
			}

			if (Power_spectrum && Lib_method)
			{
				a[id].last_hit[0] = a[id].I_input;
				a[id].last_hit[1] = a[id].m;
				a[id].last_hit[2] = a[id].h;
				a[id].last_hit[3] = a[id].n;
			}

			if (Lib_method)
				Update_Lib_way(a[id], IntOrd, 4, 4);  

			find_next_fire_time_for_single_neu(id, a[id], a_old[id], first_fire_time, t + dt);

			if (Lib_method && a_old[id].if_fired)
				a_old[id].last_fire_time += T_ref;
			update_all_neu(a, a_old, first_fire_time, t + dt);

			if (Estimate_RK4_call && a[id].t > 10)
			{
				for (int i = 0; i < N; i++)
				{
					if (Connect_Matrix[id][i] == 1)
					{
						count_num++;
						if (a[id].t - a[i].last_fire_time <= T_ref && id != i)
							syn_num++;
					}
				}
			}

		}
		else
		{
			for (int i = 0; i < N; i++)
			{
				Exchange(a[i], a_old[i]);
			}
			break;
		}
	}
}

void Generate_Poisson_times(struct neuron &a, double t, double dt)
{
	int k = a.Poisson_input_num;
	double t1;
	if (k == -1)
		t1 = -log(1 - Random(a.seed)) / a.Nu; 
	else
		t1 = a.Poisson_input_time[k];      
	k = 0;
	a.Poisson_input_time[k] = t1;
	while (t1 < t + dt)
	{
		k++;
		t1 += -log(1 - Random(a.seed)) / a.Nu;    // the last one is larger than t+dt

		if (k + 1 > int(T_Step_Large*a.Nu * 2) + 5)
			a.Poisson_input_time = (double *)realloc(a.Poisson_input_time, (k + 15) * sizeof(double));
		a.Poisson_input_time[k] = t1;
	}
	a.Poisson_input_num = k;
}

void Record_Power_spectrum(double t)
{
	fwrite(&t,sizeof(double),1,FP_FFTW);
	for (int i = 0; i < N; i++)
	{
		if (Lib_method && t - neu[i].last_fire_time < T_ref)
		{
			double output;
			Find_v_in_lib_method(neu[i], IntOrd,4,output,t-neu[i].last_fire_time);
			fwrite(&output, sizeof(double), 1, FP_FFTW);
		}
		else
			fwrite(&neu[i].v, sizeof(double), 1, FP_FFTW);
	}
}


void Run_model()
{
	double t = 0, tt = 0, tt_fftw = 0, t_lib = 0, t_fp = 0;
	double t_test = 0, s = -100, ss = neu[0].v;

	while (t < T_Max)
	{

		for (int i = 0; i < N; i++)
			Generate_Poisson_times(neu[i], t, T_Step_Large);

		evolve_model_with_initial_timestep(neu, neu_old, t, T_Step_Large);
		t += T_Step_Large;

		if (record_data[1] && t - tt >= 0.03125 && t<=2e4)
		{
			tt = t;
			fwrite(&t, sizeof(double), 1, FP1);
			for (int i = 0; i < N; i++)
			{
				double s = NAN;
				if(Lib_method && t-neu[i].last_fire_time<=T_ref)
					fwrite(&s, sizeof(double), 1, FP1);
				else
					fwrite(&neu[i].v, sizeof(double), 1, FP1);
			}
		}
		if (Power_spectrum && t <= 1e4 && t - tt_fftw >= 0.5)
		{
			tt_fftw = t;
			Record_Power_spectrum(t);
		}

		if (RecordFP && t - t_fp >= 10)				//// fire pattern 0100101010
		{
			t_fp = t;
			record_fire_pattern(t);
		}
	}
}
