void Read_parameters(long &seed, long &seed1)
{
	FILE *fp;
	fp = fopen(input_parameters, "r");

	if(fp == NULL)
	{
		printf("Error ! Cann't read parameters\n");
		getchar();// system("pause");
		exit(1);
	}

	char ch[100];

	fscanf(fp,"%s%d%s%d",ch,&NE,ch,&NI);
	N = NE+NI;	
	seed = 11, seed1 = 11;
	fscanf(fp,"%s%lf%s%lf",ch,&T_Max,ch,&T_Step_Large); 
	fscanf(fp, "%s%lf%lf%lf%lf", ch, &S[0], &S[1], &S[2], &S[3]); 
	fscanf(fp,"%s%d",ch,&I_CONST);
	fscanf(fp,"%s%lf%s%lf",ch,&Nu,ch,&f);


	fscanf(fp,"%s%lf",ch, &P_c);
	fscanf(fp, "%s%d", ch, &random_S);
	if (random_S > 4 || random_S < 0)
	{
		printf("Error pm.pS=%d\n", random_S);
		getchar();
		exit(0);
	}
	while (fgetc(fp) != '\n');

	fscanf(fp, "%s%d", ch, &random_Nu);
	if (random_Nu > 4 || random_Nu < 0)
	{
		printf("Error pm.pNu=%d\n", random_Nu);
		getchar();
		exit(0);
	}
	while (fgetc(fp) != '\n');

	fscanf(fp, "%s%d", ch, &method);
	while (fgetc(fp) != '\n');

	fscanf(fp, "%s%d", ch, &Lyapunov); 
	fscanf(fp, "%s%d%d",ch, &record_data[0], &record_data[1]);
	fscanf(fp,"%s%s%s%s",ch,file,ch,file1);  
	fclose(fp);

	if (method < 0 || method>3)
	{
		printf("Error! method=%d\n", method);
		exit(0);
	}

	if (method == 0 || method == 2 || method == 3)
		T_Step_Small = T_Step_Large;
	else
		T_Step_Small = 1.0/32;

	Power_spectrum = 0;
	Estimate_RK4_call = 0; 
	FP_fire_pattern = 0;
}


void out_put_filename()
{
	Regular_method = 0, Lib_method = 0, Adaptive_method = 0, ETDRK_method = 0;

	if (method == 0)
		Regular_method = 1;
	else if (method == 1)
		Adaptive_method = 1;
	else if (method == 2)
		Lib_method = 1;
	else if (method == 3)
		ETDRK_method = 1;
	else
	{
		printf("Error! method=%d\n", method);
		exit(0);
	}

	if (method == 0 || method == 2 || method == 3)
		T_Step_Small = T_Step_Large;

	char str[200] = "", c[10], str1[200];

	if (Lyapunov)
	{
		record_data[0] = 0;
		record_data[1] = 0;
	}
	double n, m;
	n = log(T_Step_Large) / log(0.5);
	m = log(T_Step_Small) / log(0.5);

	strcpy(str, "NE="), sprintf(c, "%d", NE), strcat(str, c);
	strcat(str, "NI="), sprintf(c, "%d", NI), strcat(str, c), strcat(str, "_");

	if (random_S == 1)
		strcat(str, "U-");
	else if (random_S == 2)
		strcat(str, "G-");
	else if (random_S == 3)
		strcat(str, "E-");
	else if (random_S == 4)
		strcat(str, "LN-");


	if(Lib_method)
	{
		strcat(str, "RK4_");
		strcat(str, "lib_"), strcat(str, "t="), sprintf(c, "%0.3f", n), strcat(str, c);
	}
	else if (ETDRK_method)
	{
		strcat(str, "ETD4RK_");
		strcat(str, "t="), sprintf(c, "%0.3f", n), strcat(str, c);
	}
	else if(Adaptive_method)
	{
		strcat(str, "RK4_");
		strcat(str, "t_l="), sprintf(c, "%0.3f", n), strcat(str, c);
		strcat(str, "t_s="), sprintf(c, "%0.3f", m), strcat(str, c);
	}
	else
	{
		strcat(str, "RK4_"),strcat(str, "t="), sprintf(c, "%0.3f", n), strcat(str, c);
	}



	strcat(str, "p="), sprintf(c, "%0.3f", P_c), strcat(str, c);
	strcat(str, "SEE="), sprintf(c, "%0.3f", S[0]), strcat(str, c);  
	strcat(str, "SIE="), sprintf(c, "%0.3f", S[1]), strcat(str, c);
	strcat(str, "SEI="), sprintf(c, "%0.3f", S[2]), strcat(str, c);
	strcat(str, "SII="), sprintf(c, "%0.3f", S[3]), strcat(str, c);
	strcat(str, "f="), sprintf(c, "%0.3f", f), strcat(str, c);
	strcat(str, "u="), sprintf(c, "%0.3f", Nu), strcat(str, c);

	if (random_Nu == 1)
		strcat(str, "U");
	else if (random_Nu == 2)
		strcat(str, "G");
	else if (random_Nu == 3)
		strcat(str, "E");
	else if (random_Nu == 4)
		strcat(str, "LN");

	printf("log2(t)=%0.2f, T_Max=%0.2ems ", -n, T_Max);

	if (record_data[0])
	{
		strcpy(str1, file), strcat(str1, str), strcat(str1, "_spikes.dat");
		FP = fopen(str1, "wb");
	}
	if (record_data[1])
	{
		strcpy(str1, file), strcat(str1, str), strcat(str1, "_voltage.dat");
		FP1 = fopen(str1, "wb");
	}

	if (Power_spectrum)				
	{
		char ch[200];
		strcpy(ch, file), strcat(ch, "fftw_"), strcat(ch, str), strcat(ch, ".dat");
		FP_FFTW = fopen(ch, "wb");
	}

	if (RecordFP)
	{
		if (N < 100)
		{
			printf("N=%d should >= 100 to record FP!\n", N);
			exit(0);
		}
		char ch[200];
		strcpy(ch, file), strcat(ch, "FP_"), strcat(ch, str), strcat(ch, ".dat");
		FP_fire_pattern = fopen(ch, "wb");
	}

	//if (record_data[0] || record_data[1] || Power_spectrum || RecordFP)
	//{
	//	printf("file:N=%d, %s\n", N, str);
	//}
	
}