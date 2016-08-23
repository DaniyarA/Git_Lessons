#include<iostream>
#include<conio.h>
#include<stdio.h>
#include<fstream>
#include<cmath>
using namespace std;

int v1[10][8]={
	//H2 O2 H2O  OH H  HO2  O  N2
	{0,  1,  0,  0, 1,  0,  0, 0},//1
	{1,  0,  0,  0, 0,  0,  1, 0},//2
	{1,  0,  0,  1, 0,  0,  0, 0},//3
	{0,  0,  1,  0, 0,  0,  1, 0},//4
	{1,  0,  0,  0, 0,  0,  0, 1},//5
	{0,  0,  0,  0, 0,  0,  2, 1},//6
	{0,  0,  0,  0, 1,  0,  0, 1},//7
};
int v2[10][8]={
   //H2  O2  H2O OH  H  HO2  O  N2
	{0,  0,  0,  1,  0,  0,  1, 0},//1
	{0,  0,  0,  1,  1,  0,  0, 0},//2
	{0,  0,  1,  0,  1,  0,  0, 0},//3
	{0,  0,  0,  2,  0,  0,  0, 0},//4
	{0,  0,  0,  0,  2,  0,  0, 1},//5
	{0,  1,  0,  0,  0,  0,  0, 1},//6
	{0,  0,  0,  0,  0,  1,  0, 1},//7
};

double dt=1e-31, eps=0.00005;
int k=0;
int k_max=4000;
double R = 8.3144598; // Дж/(моль·K)

double T = 1800.0;
double h2[500000], o2[500000], h2o[500000], oh[500000], h[500000], ho2[500000], o[500000], n2[500000];
double k_f[10], k_r[10], k_p[10];
double c[10]; // consentrasia j-go elementa


//R.S. Barlow, N.S.A. Smith  from article
				//=>			=>					=>				=>					<=				=>				<=
double A[7] = {  2e+14,			5.06e+04,			1e+08,			1.47e+10,			1.8e+18,		1.89e+15,		3.19e+18}; 
double N[7] = {  0.0,			2.67,				1.6,			1.14,				-1.0,			0.0,			-0.8	};
double Ea[7] ={  16800.0,		6286.0,				3298,			16991,				0.00,			-1788.0,		46699.3	}; //kaloriya
//REACTIONS		 H+O2<=>O+OH	O+H2<=>H+OH			H2+OH<=>H2O+H	O+H2O<=>OH+OH		H2+M<=>H+H+M	O+O+M<=>O2+M	H+O2+M<=>HO2+M

double f(int jj)
{
	double k1=0,k2=0,k3=0,k4=0,mult1,mult2;
	for(int i=0;i<10;i++)
	{
		mult1=1; mult2=1;
		for(int j=0;j<8;j++)
		{
			mult1=mult1*pow(c[j],v1[i][j]);
		}
		for(int j=0;j<8;j++)
		{
			mult2=mult2*pow(c[j],v2[i][j]);
		}
		k1=k1+(v1[i][jj]-v2[i][jj])*(k_r[i]*mult2-k_f[i]*mult1);
	}
	for(int i=0;i<10;i++)
	{
		mult1=1; mult2=1;
		for(int j=0;j<8;j++)
		{
			if(j==jj)c[j]=c[j]+dt*0.5*k1;
			mult1=mult1*pow(c[j],v1[i][j]);
		}
		for(int j=0;j<8;j++)
		{
			mult2=mult2*pow(c[j],v2[i][j]);
		}
		k2=k2+(v1[i][jj]-v2[i][jj])*(k_r[i]*mult2-k_f[i]*mult1);
	}
	for(int i=0;i<10;i++)
	{
		mult1=1; mult2=1;
		for(int j=0;j<8;j++)
		{
			if(j==jj)c[j]=c[j]+dt*0.5*k2;
			mult1=mult1*pow(c[j],v1[i][j]);
		}
		for(int j=0;j<8;j++)
		{
			mult2=mult2*pow(c[j],v2[i][j]);
		}
		k3=k3+(v1[i][jj]-v2[i][jj])*(k_r[i]*mult2-k_f[i]*mult1);
	}
	for(int i=0;i<10;i++)
	{
		mult1=1; mult2=1;
		for(int j=0;j<8;j++)
		{
			if(j==jj)c[j]=c[j]+dt*k3;
			mult1=mult1*pow(c[j],v1[i][j]);
		}
		for(int j=0;j<8;j++)
		{
			mult2=mult2*pow(c[j],v2[i][j]);
		}
		k4=k4+(v1[i][jj]-v2[i][jj])*(k_r[i]*mult2-k_f[i]*mult1);
	}
	return (c[jj]+dt*(k1+2*k2+2*k3+k4)*1.0/6);
}

int main()
{
	setlocale(LC_CTYPE, "Rus");
	ofstream outfile("res.dat", ios::out);

	//H2 O2 H2O OH H HO2 O N2 consitrations
	c[0]=0.5; c[1]=0.5; c[2]=0; c[3]=0; c[4]=0; c[5]=0; c[6]=0; c[7]=0.1;

	outfile<<"time\tH2\t\O2\t\H2O\t\OH\t\H\t\HO2\t\O\t\N2\n";

	double g[7][8], dG[7], sumG1, sumG2;
	//Dannye iz mechanisma
	double koef_a[7][8] = { // fiksirovannye koef a1,a2,..,a7 dlia kajdogo elementa 
		//H2 //O2 H2O //OH //H HO2 //O //N2
		{3.45044,		3.2129,			3.386842,		4.1253,			2.5,		3.0,		2.9464,			3.5310 },
		{-3.11E-05,		0.00112748,		0.003474,		-0.003225,		0.0,		-0.002,		-0.00163,		-0.0001236 },
		{3.30475E-07,	-5.756150E-07,	-6.3546e-06,	6.5276e-06,		0.0,		3e-06,		2.42103e-06,	-5.0299e-07 },
		{-9.2175E-11,	1.31387E-09,	6.9685e-09,		-5.798e-09,		0.0,		-3e-09,		-1.6028e-09,	2.4353e-09 },
		{7.77673E-15,	-8.7685E-13,	-2.50658e-12,	2.0623e-12,		0.0,		3e-13,		3.89096e-13,	-1.4088e-12 },
		{-1030.1098,	-1005.249,		-30208.11,		3346.30913,		25471.627,	30000,		29147.6445,		-1046.976 },
		{-3.959938,		6.03473,		2.59023,		-0.69043296,	-0.46011,	3.000,		2.963994,		2.96747 }
	};

	do {

		k_f[0] = (3.0e+14) * exp((-8.81)/( T));		/*H+O2 = O+OH  */	k_r[0] = 2.48e+13*exp(-0.66/T);
		k_f[1] = (3.0e+14) * exp((-4.03)/( T));		/*O+H2=H+OH    */	k_r[1] = 1.3e+14*exp(-2.49/T);
		k_f[2] = (3.0e+14) * exp((-3.02)/( T));		/*H2+OH=H2O+H  */	k_r[2] = 1.33e+15*exp(-10.95/T);
		k_r[3] = (3.0e+14) * exp((-3.02/T)/( T));	/*O+H2O=OH+OH  */	k_f[3] = 3.12e+15*exp(-12.51/T);
		k_f[4] = (1.35e+17) * exp((-54.0)/( T)/T);	/*H2+M=H+H+M   */	k_r[4] = 1.0e+16;
		k_r[5] = (5.8e+16) * exp((-60.6/T)/( T));	/*O+O+M=O2+M   */	k_f[5] = 6.0e+14;
		k_r[6] = (8.0e+16) * exp((-52.0/T)/( T));	/*O+H+M=OH+M   */	k_f[6] = 1.0e+16;
		k_r[7] = (9.66e+17) * exp((-54.0/T)/( T));	/*H+OH+M=H2O+M */	k_f[7] = 1.0E+17;
		k_r[8] = (2.4e+15) * exp((-23.1)/( T));		/*H+O2+M=HO2+M */	k_f[8] = 1.59e+15 * exp(0.203/T);
		//k_f[9] = (2.5e+14) * exp((-0.95)/( T));		/*HO2+H=OH+OH  */	k_r[9] = 1.20e+13 * exp(-20.2/T);
		/*
		printf("************** Выходные данные *********************\n");
		printf("k_fi\t\t\t\t");
		printf("dGi\t\t\t\t");
		printf("k_pi\t\t\t\t");
		printf("k_ri\n");

		for(int i = 0; i<7; i++)
		{
			//******************** k_fi скорость прямой реакции **************************************

			k_f[i] = A[i] * pow(T, N[i]) * exp((-1.0 * 4.184 * Ea[i])/(R * T)); // dlia par-v iz Workbanch
			//k_f[i] = ((A[i])/pow(10, 6)) * pow(T, N[i]) * exp(-1.0*Ea[i]/(R * T)); // dlia pr iz Varnarts
			printf("%.15E\t\t", k_f[i]);

			//******************** g[i][j] *************
			for(int j = 0; j<8;j++){
				g[i][j]=R*(koef_a[0][j]*(T - T*log(T)) + 0.5*koef_a[1][j]*pow(T,2) + koef_a[2][j]*(1.0/6.0)*pow(T,3) + koef_a[3][j]*(1.0/12.0)*pow(T,4) + koef_a[4][j]*(1.0/20.0)*pow(T, 5) + koef_a[5][j] - koef_a[6][j]*T);
			}

			//******************** dG[i] ***************
			sumG1 = 0.0;
			for(int j = 0; j<8;j++){
				sumG1 = sumG1 + v1[i][j]*g[i][j];
			}

			sumG2 = 0.0;
			for(int j = 0; j<8;j++){
				sumG2 = sumG2 + v2[i][j]*g[i][j];
			}

			dG[i] = sumG2 - sumG1;
			printf("%.15E\t\t", dG[i]);

			//******************** k_pi скорость прямой реакции **************************************

			k_p[i] = exp(-1.0*dG[i]/(R*T));
			printf("%.15E\t\t", k_p[i]);

			//******************** k_ri скорость прямой реакции **************************************

			k_r[i]=k_f[i]/k_p[i];
			printf("%.15E\t\t", k_r[i]);

			printf("\n");
		}
		//******************* runge-kutte ****************
		*/
		for(int j = 0; j<8;j++)
		{
			c[j]=f(j); //runge-kutte v func
		}

		//if(k%100==0)
		//{
		printf("%.5E\t", c[0]); printf("%.7E\t", c[1]);
		printf("%.5E\t", c[2]); printf("%.5E\t", c[3]);
		printf("%.5E\n\n", c[4]);
		outfile<<k<<"\t"<<c[0]<<"\t"<<c[1]<<"\t"<<c[2]<<"\t"<<c[3]<<"\t"<<c[4]<<"\t"<<c[5]<<"\t"<<c[6]<<"\t"<<c[7]<<endl;
		//}

		k++;
	}while( k<1000);



	printf("\n");

	for(int i = 0; i<10; i++)
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

