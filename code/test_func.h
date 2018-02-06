int neu_id[10] = { 88,21,29,11,23,81,80,24,44,85 };


///// state pattern, 010100010101
void record_fire_pattern(double t)
{
	double s[10];
	for (int i = 0; i < 10; i++)
		if (t - neu[neu_id[i]].last_fire_time < 10)
			s[i] = 1;
		else
			s[i] = 0;	
	fwrite(s, sizeof(double), 10, FP_fire_pattern);
}



