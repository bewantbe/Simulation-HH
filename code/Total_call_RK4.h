double new_total_call_RK4(double R)
{
	double num, dt, T, dts;
	double p, p1;

	dt = T_Step_Large;
	dts = T_Step_Small;
	p = T_ref*R;
	T = T_Max;

	if (count_num == 0)
		p1 = 0;
	else
		p1 = syn_num / count_num;

	if (Lib_method)
	{
		num = (1 - p)*(T / dt + T*Nu) + T*(2 + Nu*dt)*(R + R*P_c*(N - 1)*(1 - p1));
	}
	else if (Adaptive_method)
	{
		num = T / dt + T*Nu + T*(2 + Nu*dt)*(R + R*P_c*(N - 1)) + T*R*(T_ref / dts + 1) + T*R*P_c*(N - 1)*p1*dt / 2 / dts;
	}
	else
	{
		num = T / dt + T*Nu + T*(2 + Nu*dt)*(R + R*P_c*(N - 1)); 
	}
	return num*N;
}