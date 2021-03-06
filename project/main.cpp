#include<iostream>
#include<conio.h>
#include<stdio.h>
#include<cmath>
using namespace std;
int main()
{
	setlocale(LC_CTYPE, "Rus");
	/*
	// iz knigi Varnartsa
	double A[7] = {2E+14, 5.06E+04, 1.0E+08, 1.5E+09, 1.8E+18, 2.2E+22,2.5E+13}; 
	double N[7] = {0, 2.670, 1.6, 1.14, -1.0, -2.0, 0.0};
	double Ea[7] = {70.3E+03, 26.3E+03, 13.8E+03, 0.42E+03, 0.0E+03, 0.0E+03, 2.9E+03}; //joul
	*/

	//iz workbancha
	double A[7] = {5.88993E-09, 8.43554E-20, 3.58676E-16, 4.9318E-18, 7.60029E-05, 1.04781E-25, 2.75649E-11}; 
	double N[7] = {-0.406, 2.670, 1.510, 2.02, -1.400, -2.0, 0.600};
	double Ea[7] = {16.599E+03, 6.29E+03, 3.43E+03, 13.40E+03, 104.38E+03, 0.0E+03, 0.823E+03}; //kallorii

	double R = 8.3144598; // ��/(�����K)
	double k_f[10], k_r[10], k_p[10];

	int v1[7][8]={
		//H2 O2 H2O OH H HO2 O N2
		{0, 1, 0, 0, 1, 0, 0, 0},
		{1, 0, 0, 0, 0, 0, 1, 0},
		{1, 0, 0, 1, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0, 1, 0},
		{1, 0, 0, 0, 0, 0, 0, 1},
		{0, 0, 0, 1, 1, 0, 0, 1},
		{0, 0, 0, 0, 1, 1, 0, 0}
	};

	int v2[7][8]={
		//H2 O2 H2O OH H HO2 O N2
		{0, 0, 0, 1, 0, 0, 1, 0},
		{0, 0, 0, 1, 1, 0, 0, 0},
		{0, 0, 1, 0, 1, 0, 0, 0},
		{0, 0, 0, 2, 0, 0, 0, 0},
		{0, 0, 0, 0, 2, 0, 0, 1},
		{0, 0, 1, 0, 0, 0, 0, 1},
		{1, 1, 0, 0, 0, 0, 0, 0}
	};

	printf("************** ������� ������ *********************\n");
	// vyvod znachenii A, b, Ea, V1,V2
	printf("A_i:\n");
	for(int i=0; i<7; i++)
		printf("\t%.5E\n", A[i]);
	printf("\n");

	printf("N_i:\n");
	for(int i=0; i<7; i++)
		printf("\t%.5E\n", N[i]);
	printf("\n");

	printf("Ea_i:\n");
	for(int i=0; i<7; i++)
		printf("\t%.5E\n", Ea[i]);
	printf("\n");

	printf("v1:\n");
	printf("\tH2 O2 H2O OH H HO2 O N2\n");
	for(int i = 0; i <7 ;i++){
		printf("\t");
		for(int j=0;j<8;j++){
			printf("%d ", v1[i][j]);
		}
		printf("\n");
	}
	printf("\nv2:\n");
	printf("\tH2 O2 H2O OH H HO2 O N2\n");
	for(int i = 0; i <7 ;i++){
		printf("\t");
		for(int j=0; j<8;j++){
			printf("%d ", v2[i][j]);}
		printf("\n");
	}

	double T = 2000.0;

	printf("\n");
	printf("\nT = %f\n\n", T);

	double g[7][8], dG[7], sumG1, sumG2;
	//Dannye iz mechanisma
	double koef_a[7][8] = { // fiksirovannye koef a1,a2,..,a7 dlia kajdogo elementa 
		//H2 //O2 H2O //OH //H HO2 //O //N2
		{3.45044, 3.2129, 3.386842, 4.1253, 2.5, 3.0, 2.9464, 3.5310 },
		{-3.11E-05, 0.00112748, 0.003474, -0.003225, 0.0, -0.002, -0.00163, -0.0001236 },
		{3.30475E-07, -5.756150E-07, -6.3546e-06, 6.5276e-06, 0.0, 3e-06, 2.42103e-06, -5.0299e-07 },
		{-9.2175E-11, 1.31387E-09, 6.9685e-09, -5.798e-09, 0.0, -3e-09, -1.6028e-09, 2.4353e-09 },
		{7.77673E-15, -8.7685E-13, -2.50658e-12, 2.0623e-12, 0.0, 3e-13, 3.89096e-13, -1.4088e-12 },
		{-1030.1098, -1005.249, -30208.11, 3346.30913, 25471.627, 30000, 29147.6445, -1046.976 },
		{-3.959938, 6.03473, 2.59023, -0.69043296, -0.46011, 3.000, 2.963994, 2.96747 }
	};

	/*
	//Dannye iz BD Workbench
	double koef_a[7][8] = { // fiksirovannye koef a1,a2,..,a7 dlia kajdogo elementa 
	//H2 //O2 H2O //OH //H HO2 //O //N2
	{3.45044, 3.1340E+000, 3.9043E+000, 3.8539, 2.5, 3.298, 2.7961, 3.2658 },
	{-3.11E-05, 1.5051E-003, 1.2084E-005, -1.3381E-003, 0.0, 2.8775E-03, -7.0062E-004, 8.5495E-004 },
	{3.30475E-07, -5.8959E-007, 1.7109E-006, 1.7642E-006, 0.0, 1.1055E-06, 6.0328E-007, -2.0802E-007 },
	{-9.2175E-11, 1.1509E-010, -7.4481E-010, -6.7959E-010, 0.0, -2.3583E-9, -2.1871E-010, 2.0869E-011 },
	{7.77673E-15, -8.1640E-015, 8.8070E-014, 8.7096E-014, 0.0, 8.0648E-13, 2.8260E-014, -6.2642E-016 },
	{-1.0301E+003, -9.9633E+002, -3.0263E+004, 3.6287E+003, 2.5474E+4, 4.8625E+1, 2.9162E+004, -1.0099E+003 },
	{-3.9599E+000, 6.3803E+000, 3.8002E-001, 4.5376E-001, -0.45995, 7.8628, 3.6105E+000, 4.1786E+000 }
	};
	*/
	printf("************** �������� ������ *********************\n");
	printf("k_fi\t\t\t\t");
	printf("dGi\t\t\t\t");
	printf("k_pi\t\t\t\t");
	printf("k_ri\n");
	for(int i = 0; i<7; i++)
	{
		//******************** k_fi �������� ������ ������� **************************************

		k_f[i] = (1.0/A[i]) * pow(T, N[i]) * exp((-1.0 * 4.1840 * Ea[i])/(R * T)); // dlia par-v iz Workbanch
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

		//******************** k_pi �������� ������ ������� **************************************

		k_p[i] = exp(-1.0*dG[i]/(R*T));
		printf("%.15E\t\t", k_p[i]);

		//******************** k_ri �������� ������ ������� **************************************

		k_r[i]=k_f[i]/k_p[i];
		printf("%.15E\t\t", k_r[i]);

		printf("\n");
	}

	//******************* 
	printf("\n");
	system("pause");
	return 0;
}

