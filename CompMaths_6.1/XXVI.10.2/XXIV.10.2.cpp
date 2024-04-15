#include <iostream>
#include <cmath>
#include <fstream>

const double M_PI = 3.14159265359;

double L = 20;
double current_number = 0.6;
double h = 0.5;
int NX = L/h;
int T = 18;
double t = current_number * h;
int T_counts = T / t;

void Solver()
{

	double** SolMat = new double* [T_counts];
	for (int i = 0; i < T_counts; ++i)
	{
		SolMat[i] = new double[NX];
	}

	for (int i = 0; i < NX; ++i)
	{
		SolMat[0][i] = sin(2 * M_PI * h * i / L);
	}

	SolMat[0][NX-1] = SolMat[0][0];

	for (int i = 0; i < T_counts - 1; ++i)
	{
		for (int j = NX-1; j >= 0; --j)
		{
			if (j == 0)
				SolMat[i + 1][0] = -current_number * (SolMat[i][0] - SolMat[i][NX - 1]) + SolMat[i][0];
			else
				SolMat[i + 1][j] = -current_number * (SolMat[i][j] - SolMat[i][j - 1]) + SolMat[i][j];
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
	for (int i = 0; i < NX - 1; ++i)
	{
		SolMat[0][i] = sin(2 * M_PI * h * i / L);
	}
	SolMat[0][NX - 1] = SolMat[0][0];

	for (int i = 0; i < T_counts - 1; ++i)
	{
		for (int j = NX; j >= 0; --j)
		{
			if (j == 0)
				SolMat[i + 1][0] = -current_number * current_number/2*
				(SolMat[i][1]- SolMat[i][0] + SolMat[i][NX - 1]) -
				current_number/2*(SolMat[i][1]);
			else
				SolMat[i + 1][j] = -current_number * (SolMat[i][j] - SolMat[i][j - 1]) + SolMat[i][j];
		}
	}
}

int main()
{

	Solver();
	return 0;
}