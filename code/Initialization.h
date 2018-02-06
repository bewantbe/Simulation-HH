double Decide_S(int i, int j, long &seed, long &seed1)  // i-->j
{
	if (i < NE && j < NE)
	{
		double s = S[0]; 
		if (random_S == 0)
			return s;
		else if (random_S == 1)						// uniform [0, 2s]
			return Random(seed) * 2 * s;
		else if (random_S == 2)						// gauss N(s,(s/4)^2)
			return abs(sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*s / 4 + s);
		else if (random_S == 3)                     // Exponential E(s)
			return -log(1 - Random(seed)) * s;
		else										// Log normal mu = 1.25*log(s), sigma^2 = -log(s)/2
		{
			double a = log(s);
			double b = sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*sqrt(-a / 2) + 1.25*a;
			return exp(b);
		}
	}
	else if (i < NE && j >= NE)
	{
		double s = S[1]; 
		if (random_S == 0)
			return s;
		else if (random_S == 1)
			return Random(seed) * 2 * s;
		else if (random_S == 2)
			return abs(sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*s / 4 + s);
		else if (random_S == 3)                    
			return -log(1 - Random(seed)) * s;
		else										
		{
			double a = log(s);
			double b = sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*sqrt(-a / 2) + 1.25*a;
			return exp(b);
		}
	}
	else if (i >= NE && j < NE)
	{
		double s = S[2];		
		if (random_S == 0)
			return s;
		else if (random_S == 1)
			return Random(seed) * 2 * s;
		else if (random_S == 2)
			return abs(sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*s / 4 + s);
		else if (random_S == 3)
			return -log(1 - Random(seed)) * s;
		else
		{
			double a = log(s);
			double b = sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*sqrt(-a / 2) + 1.25*a;
			return exp(b);
		}
	}
	else
	{
		double s = S[3];		
		if (random_S == 0)
			return s;
		else if (random_S == 1)
			return Random(seed) * 2 * s;
		else if (random_S == 2)
			return abs(sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*s / 4 + s);
		else if (random_S == 3)
			return -log(1 - Random(seed)) * s;
		else
		{
			double a = log(s);
			double b = sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*sqrt(-a / 2) + 1.25*a;
			return exp(b);
		}
	}
}

double Decide_Nu(int i, long &seed, long &seed1)  // i-->j
{
	double s = Nu;
	if (random_Nu == 0)
		return s;
	else if (random_Nu == 1)					 // uniform [0, 2s]
		return Random(seed) * 2 * s;
	else if (random_Nu == 2)					 // gauss N(s,(s/4)^2)
		return abs(sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*s / 4 + s);
	else if (random_Nu == 3)                     // Exponential E(s)
		return -log(1 - Random(seed)) * s;
	else										 // Log normal, log(X)~N(mu,sigma^2)
	{
		double a = log(s);
		double b = sqrt(-2 * log(Random(seed)))*cos(2 * PI*Random(seed1))*sqrt(-a / 2) + 1.25*a;
		return exp(b);
	}	
}

void Create_connect_matrix(long& seed)
{
	Connect_Matrix = new double *[N];
	for (int i = 0; i < N; i++)
		Connect_Matrix[i] = new double[N];

	CS = new double *[N];
	for (int i = 0; i < N; i++)
		CS[i] = new double[N];
	long seed1 = 125;	      

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			Connect_Matrix[i][j] = 0;
			CS[i][j] = 0;
		}


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j && Random(seed) < P_c)
			{
				Connect_Matrix[i][j] = 1;
				CS[i][j] = Decide_S(i, j, seed, seed1);
			}
		}
	}


	if (record_data[0] || record_data[1])
	{
		FILE *fp;
		char str[200], ch[10], c[10];

		strcpy(str, file);
		strcat(str, "NE="), sprintf(c, "%d", NE), strcat(str, c);
		strcat(str, "NI="), sprintf(c, "%d", NI), strcat(str, c), strcat(str, "_");
		strcat(str, "connect_matrix-p=");
		sprintf(ch, "%0.3f", P_c), strcat(str, ch);

		if (random_S == 1)
			strcat(str, "-U");
		else if (random_S == 2)
			strcat(str, "-G");
		else if (random_S == 3)
			strcat(str, "-E");
		else if (random_S == 4)
			strcat(str, "-LN");
		strcat(str, ".dat");

		if ((fp = fopen(str, "rb")) == NULL)
		{
			fp = fopen(str, "wb");
			for (int i= 0; i < N; i++)
				fwrite(Connect_Matrix[i], sizeof(double), N, fp);

			if (random_S != 0)
				for (int i = 0; i < N; i++)
					fwrite(CS[i], sizeof(double), N, fp);
			fclose(fp);
		}

	}
}

void Find_unique_resolution_equalspace(double *a, int n, int k,int &Resolution)  
{
	double *b = new double[n];
	double s = 0,s1;
	for (int i = 0; i < n; i++)
		b[i] = a[i];
	sort(b, b + n);
	unique(b, b + n);
	Resolution = 1;
	for (int i = 1; i < n; i++)
	{
		if (b[i] <= b[i - 1])
			break;
		else
			Resolution++;
	}

	Lib_unique[k] = new double[Resolution];
	for (int i = 0; i < Resolution; i++)
	{
		Lib_unique[k][i] = b[i];
		s += b[i];
	}

	s1 = (b[0] + b[Resolution - 1])*Resolution / 2.0;
	int equal_space = abs(s1 - s) < Epsilon ? 1 : 0;
	delete[] b;
	if(!equal_space)
	{
		printf("Error! Not eaual sapced sampling: %d \n ",k);
		getchar();
		//system("pause");
		exit(1);
	}
}

void Find_unique_resolution_equalspace_new(double *a, int n, int k, int &Resolution)
{
	double min, max, step;

	min = a[0], max = a[0];
	for (int i = 1; i < n; i++)
	{
		if (a[i] > max)
			max = a[i];
		if (a[i] < min)
			min = a[i];
	}

	double sec_min = max;
	for (int i = 1; i < n; i++)
	{
		if (a[i]<sec_min && a[i]>min)
			sec_min = a[i];
	}
	step = sec_min - min;

	Resolution = int((max - min) / step + 1.001);

	Lib_unique[k] = new double[Resolution];
	for (int i = 0; i < Resolution; i++)
		Lib_unique[k][i] = min+i*step;
}

void Initialize_library()
{
	char str[200];
	int length = 0, length_old = 1.2e5, n = 8;
	clock_t t0, t1;

	strcpy(str, file1), strcat(str, "Lib.dat");

	FILE *fp = fopen(str, "rb");
	if (fp == NULL)
	{
		printf("Error in Initialization()! :: Cann't open Lib.dat! \n");
		getchar();// system("pause");
		exit(1);
	}

	Lib_data = new double*[n];
	for (int i = 0; i < n; i++)
		Lib_data[i] = new double[length_old];
	Lib_unique = new double*[4]; //I_input,m,h,n


	double *s = new double[n];
	int read = 1;
	while (!feof(fp))
	{
		if (length + 1 > length_old)
		{
			for (int i = 0; i < n; i++)
				Lib_data[i] = (double*)realloc(Lib_data[i], (length + 1e5) * sizeof(double));
			length_old = length + 1e5;
		}


		if (fread(s, sizeof(double), n, fp) != n) // dat file
		{
			read = 0;    // End of file (EOF)
			break;
		}
		if (read)
		{
			for (int i = 0; i < n; i++)
				Lib_data[i][length] = s[i];
			length++;
		}
	}
	delete[] s;
	fclose(fp);

	Lib_length = length;
//	printf("Lib_length = %d\n", Lib_length);

	Lib_resolution = new int[4];		//I_input,m,h,n

	for (int i = 0; i < 4; i++)
		Find_unique_resolution_equalspace(Lib_data[i], Lib_length, i, Lib_resolution[i]);
	
}


// if power spectrum, trace for v
void Initial_library_trace()
{
	char str[200];
	int length = 0, length_old = 1.2e5, n = int(T_ref*256+1+0.1);

	strcpy(str, file1), strcat(str, "Lib_v.dat");

	FILE *fp;   //v,m,h,n

	fp = fopen(str, "rb");
	if (fp == NULL)
	{
		printf("Error! :: Cann't open Lib_v ! \n");
		getchar();// system("pause");
		exit(1);
	}


	Lib_v = new double*[n];
	for (int i = 0; i < n; i++)
		Lib_v[i] = new double[length_old];


	double *s = new double[n];

	int read = 1;
	while (!feof(fp))
	{
		if (length + 1 > length_old)
		{
			for (int i = 0; i < n; i++)
				Lib_v[i] = (double*)realloc(Lib_v[i], (length + 1e5) * sizeof(double));
			length_old = length + 1e5;
		}

		if (fread(s, sizeof(double), n, fp) != n) // dat file
		{
			read = 0;    // End of file (EOF)
			break;
		}

		if (read)
		{
			for (int i = 0; i < n; i++)
				Lib_v[i][length] = s[i];
			length++;
		}
	}


	delete[]s;
	fclose(fp);

	if (length != Lib_length)
	{
		printf("Error in length! Lib_length=%d Lib_v_length=%d\n", Lib_length, length);
		getchar();
		exit(0);
	}
}

void Initialization(long &seed0,long &seed2)
{
	if (Lib_method)
		Initialize_library();

	if (Lib_method && Power_spectrum)
		Initial_library_trace();

	Create_connect_matrix(seed0);
	neu = new struct neuron[N];
	neu_old = new struct neuron[N];


	long Seed = 11, Seed1 = 125;

	for (int i = 0; i < N; i++)
	{
		neu[i].t = 0;
		neu[i].Nu = Decide_Nu(i,Seed,Seed1);
		neu[i].v = -65;
		neu[i].dv = 0;
		neu[i].m = alpha_m(neu[i].v) / (alpha_m(neu[i].v) + beta_m(neu[i].v));
		neu[i].h = alpha_h(neu[i].v) / (alpha_h(neu[i].v) + beta_h(neu[i].v));
		neu[i].n = alpha_n(neu[i].v) / (alpha_n(neu[i].v) + beta_n(neu[i].v));
		neu[i].G_se = 0;
		neu[i].G_sse = 0;
		neu[i].G_si = 0;
		neu[i].G_ssi = 0;
		neu[i].G_f = 0;
		neu[i].G_ff = 0;
		neu[i].I_input = 0;
		neu[i].fire_num = 0;
		neu[i].last_fire_time = -1e5;
		neu[i].if_fired = 0;
		neu[i].Poisson_input_time = new double[int(T_Step_Large*Nu*2)+5];
		for (int i = 0; i < 500; i++)
			Random(seed2); 
		neu[i].seed = seed2;
		neu[i].Poisson_input_num = -1;

		neu[i].wait_strength_E = 0;
		neu[i].wait_strength_I = 0;

		neu[i].id_F = 0;
		neu[i].state = 1;
		neu_old[i].state = 0;
	}
	
}
