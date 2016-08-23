#include<iostream>
#include<conio.h>
#include<stdio.h>
#include<fstream>
#include<cmath>
using namespace std;

int v1[10][8]={
	//H2 O2 H2O OH H HO2 O N2
	{0, 1, 0, 0, 1, 0, 0, 0},//1
	{1, 0, 0, 0, 0, 0, 1, 0},//2
	{1, 0, 0, 1, 0, 0, 0, 0},//3
	{0, 0, 0, 2, 0, 0, 0, 0},//4
	{1, 0, 0, 0, 0, 0, 0, 1},//5
	{0, 0, 1, 0, 0, 0, 0, 1},//6
	{0, 0, 0, 1, 0, 0, 0, 1},//7
	{0, 1, 0, 0, 0, 0, 0, 1},//8
	{1, 1, 0, 0, 0, 0, 0, 0}//9
};
int v2[10][8]={
	//H2 O2 H2O OH H HO2 O N2
	{0, 0, 0, 1, 0, 0, 1, 0},//1
	{0, 0, 0, 1, 1, 0, 0, 0},//2
	{0, 0, 1, 0, 1, 0, 0, 0},//3
	{0, 0, 1, 0, 0, 0, 1, 0},//4
	{0, 0, 0, 0, 2, 0, 0, 1},//5
	{0, 0, 0, 1, 1, 0, 0, 1},//6
	{0, 0, 0, 0, 1, 0, 1, 1},//7
	{0, 0, 0, 0, 0, 0, 2, 1},//8
	{0, 0, 0, 2, 0, 0, 0, 0}//9
};

double dt=1e-19, eps=0.00005;
int k=0;
int k_max=4000;
double R = 8.3144598; // Дж/(моль·K)

double T = 1700.0;
double k_f[10], k_r[10], k_p[10];
double c[10]; // consentrasia j-go elementa

double RungeKuttaFunc(int jj)
{
	double k1=0, k2=0, k3=0, k4=0, mult1, mult2;
	for(int i=0; i<10; i++)
	{
		mult1 = 1; 
		mult2 = 1;
		for(int j=0; j<8; j++)
		{
			mult1 = mult1 * pow(c[j], v1[i][j]);
		}
		for(int j = 0; j<8; j++)
		{
			mult2 = mult2 * pow(c[j], v2[i][j]);
		}
		k1 = k1 + (v2[i][jj] - v1[i][jj]) * (k_f[i] * mult1 - k_r[i]*mult2);
	}
	for(int i=0; i<10; i++)
	{
		mult1 = 1; 
		mult2 = 1;
		for(int j=0; j<8; j++)
		{
			if(j==jj) 
				c[j] = c[j] + 0.5 * dt * k1;
			mult1 = mult1 * pow(c[j], v1[i][j]);
		}
		for(int j=0; j<8; j++)
		{
			mult2 = mult2 * pow(c[j], v2[i][j]);
		}
		k2 = k2 + (v2[i][jj] - v1[i][jj]) * (k_f[i] * mult1 - k_r[i] * mult2);
	}  
	for(int i=0; i<10; i++)
	{
		mult1 = 1; 
		mult2 = 1;
		for(int j=0; j<8; j++)
		{
			if(j==jj)
				c[j] = c[j] + 0.5 * dt * k2;
			mult1 = mult1 * pow(c[j], v1[i][j]);
		}
		for(int j=0; j<8; j++)
		{
			mult2 = mult2 * pow(c[j], v2[i][j]);
		}
		k3 = k3 + (v2[i][jj] - v1[i][jj]) * (k_f[i] * mult1 - k_r[i] * mult2);
	}
	for(int i=0; i<10; i++)
	{
		mult1 = 1; 
		mult2 = 1;
		for(int j=0; j<8; j++)
		{
			if(j==jj)
				c[j] = c[j] + dt * k3; 
			mult1 = mult1 * pow(c[j], v1[i][j]);
		}
		for(int j = 0; j<8; j++)
		{
			mult2 = mult2 * pow(c[j], v2[i][j]);
		}
		k4 = k4 + (v2[i][jj] - v1[i][jj]) * (k_f[i] * mult1 - k_r[i] * mult2);
	}
	return (c[jj] + dt * (k1 + 2 * k2 + 2 * k3 + k4) * (1.0/6.0));
}

int main()
{
	setlocale(LC_CTYPE, "Rus");
	ofstream outfile("resEA9reac.dat", ios::out);

	//H2 O2 H2O OH H HO2 O N2 consitrations
	//c[0]=1.0; c[1]=0.266; c[2]=0.256; c[3]=0; c[4]=0; c[5]=0; c[6]=0.2; c[7]=0.478;
	c[0]=1.0; c[1]=0.99; c[2]=0; c[3]=0; c[4]=0; c[5]=0; c[6]=0; c[7]=0.01;

	outfile<<"time\tH2\t\O2\t\H2O\t\OH\t\H\t\HO2\t\O\t\N2\n";

	do {
		k_f[0] = (3.0e+14) * exp((-8.81)/( T));			/* H+O2 = O+OH   */		k_r[0] = 2.48e+13*exp(-0.66/T);
		k_f[1] = (3.0e+14) * exp((-4.03)/( T));			/* O+H2=H+OH     */		k_r[1] = 1.3e+14*exp(-2.49/T);
		k_f[2] = (3.0e+14) * exp((-3.02)/( T));			/* H2+OH=H2O+H   */		k_r[2] = 1.33e+15*exp(-10.95/T);
		k_f[3] = (3.0e+14) * exp((-3.02/T)/( T));		/* O+H2O=OH+OH   */		k_r[3] = 3.12e+15*exp(-12.51/T);
		k_f[4] = (1.35e+17) * exp((-54.0)/( T)/T);		/* H2+M=H+H+M    */		k_r[4] = 1.0e+16;
		k_f[5] = (9.66e+17) * exp((-54.0/T)/( T));		/* H+OH+M=H2O+M  */		k_r[5] = 1.0E+17;
		k_f[6] = (8.0e+16) * exp((-52.0/T)/( T));		/* O+H+M=OH+M    */		k_r[6] = 1.0e+16;
		k_f[7] = (5.8e+16) * exp((-60.6/T)/( T));		/* O+O+M=O2+M    */		k_r[7] = 6.0e+14;
		k_f[8] = (2.72e+13) * exp((-29.2)/( T));		/* H+O2+M=HO2+M  */		k_r[8] = 8.34e+10 * exp(-1.032/T);
		// k_f[9] = (2.5e+14) * exp((-0.95)/( T)); /*HO2+H=OH+OH */ k_r[9] = 1.20e+13 * exp(-20.2/T);


		//******************* runge-kutte ****************

		if(k%100==0)
		{
		printf("%.5E\t", c[0]); 
		printf("%.7E\t", c[1]);
		printf("%.5E\t", c[2]); 
		printf("%.5E\t", c[3]);
		printf("%.5E\n", c[4]);
		printf("%.7E\t", c[5]);
		printf("%.5E\t", c[6]); 
		printf("%.5E\t", c[7]);
		outfile<<k<<"\t"<<c[0]<<"\t"<<c[1]<<"\t"<<c[2]<<"\t"<<c[3]<<"\t"<<c[4]<<"\t"<<c[5]<<"\t"<<c[6]<<"\t"<<c[7]<<endl;
		}

		for(int j = 0; j<8;j++)
		{
			c[j]=RungeKuttaFunc(j); //runge-kutte v func
		}
		k++;
	}while( k<100000);



	printf("\n");

	for(int i = 0; i<9; i++)
	{ //vyvod k_f i k_r
		printf("%.5E\t", k_f[i]);
		printf("%.5E\t", k_r[i]);
		printf("\n");
	}
	outfile.close();
	printf("\n");
	system("pause");
	return 0;
}