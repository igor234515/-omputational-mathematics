#include <iostream>
#include <cmath>
#include <fstream>

const double M_PI = 3.14159265359;

double L = 20;
double CFL = 0.6;
double h = 0.5;
int NX = L/h;
int T = 18;
double t = CFL * h;
int T_counts = T / t;

void Solver()
{

	double** SolMat = new double* [T_counts];
	for (int i = 0; i < T_counts; ++i)
	{
		SolMat[i] = new double[NX];
	}

	for (int i = 0; i <= NX-1; ++i)
	{
		SolMat[0][i] = sin(4 * M_PI * h * i / L);
	}

	SolMat[0][NX-1] = SolMat[0][0];

	

	for (int i = 0; i < T_counts - 1; ++i)
	{
		for (int j = NX-1; j >= 0; --j)
		{
			if (j == 0)
				SolMat[i + 1][0] = -CFL * (SolMat[i][0] - SolMat[i][NX - 1]) + SolMat[i][0];
			else
				SolMat[i + 1][j] = -CFL * (SolMat[i][j] - SolMat[i][j - 1]) + SolMat[i][j];
		}
	}

	std::ofstream outfile;
	outfile.open("Solution.txt", std::ios_base::out | std::ios_base::trunc);

	for (int i = 0; i < T_counts; ++i)
	{
		for (int j = 0; j < NX; ++j)
		{
			outfile << SolMat[i][j] << " ";
		}
		outfile << std::endl;
	}
	outfile.close();


	delete[] SolMat;
}

void Laks_Wendorf()
{
	double** SolMat = new double* [T_counts];
	for (int i = 0; i < T_counts; ++i)
	{
		SolMat[i] = new double[NX];
	}
	for (int i = 0; i <= NX-1; ++i)
	{
		SolMat[0][i] = sin(4 * M_PI * h * i / L);
	}
	SolMat[0][NX - 1] = SolMat[0][0];



	for (int i = 0; i < T_counts - 1; ++i)
	{
		for (int j = 1; j <= NX-2; ++j)
		{
			SolMat[i + 1][j] = CFL * CFL / 2 *
				(SolMat[i][j + 1] - 2 * SolMat[i][j] + SolMat[i][j - 1]) -
				CFL / 2 * (SolMat[i][j + 1] - SolMat[i][j - 1]) + SolMat[i][j];
		}
			SolMat[i + 1][0] = CFL * CFL / 2 *
			(SolMat[i][1] - 2 * SolMat[i][0] + SolMat[i][NX - 1]) -
			CFL / 2 * (SolMat[i][1] - SolMat[i][NX - 1]) + SolMat[i][0];

			SolMat[i + 1][NX - 1] = CFL * CFL / 2 *
			(SolMat[i][0] - 2 * SolMat[i][NX - 1] + SolMat[i][NX - 2]) -
			CFL / 2 * (SolMat[i][0] - SolMat[i][NX - 2]) + SolMat[i][NX - 1];
	}

	std::ofstream outfile2;
	outfile2.open("Laks.txt", std::ios_base::out | std::ios_base::trunc);

	for (int i = 0; i < T_counts; ++i)
	{
		for (int j = 0; j <= NX-1; ++j)
		{
			outfile2 << SolMat[i][j] << " ";
		}
		outfile2 << std::endl;
	}
		outfile2.close();
	

	delete[] SolMat;
}

int main()
{
	Solver();
	//Laks_Wendorf();
	return 0;
}