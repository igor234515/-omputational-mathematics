#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <fstream>

//general
const double L = 10.0;
const double T = 0.015;
const double γ = 5.0 / 3.0;
const int NX = 100;
double h = 2.0 * L / NX;
const double CFL = 0.01;
//left
double ρL = 13;
double PL = 1010000.0;
double uL = 0;
//right
double ρR = 1.3;
double PR = 101000.0;
double uR = 0;


double P(double ρe)
{
	return (γ - 1) * ρe;
}

double c(double ρ, double ρe)
{
	return sqrt(γ * P(ρe) / ρ);
}

double u(double ρ, double ρu)
{
	return ρu / ρ;
}


void Solver()
{

	double ρeL = PL / (γ - 1.0);
	double ρeR = PR / (γ - 1.0);

	Eigen::Matrix<double, 3, NX> w;
	Eigen::Matrix<double, 3, NX> w1;
	for (int i = 0; i < NX/2; ++i)
	{
		w(0, i) = ρL;
		w(1, i) = ρL*uL;
		w(2, i) = ρeL;
	}

	for (int i = NX/2; i < NX; ++i)
	{
		w(0, i) = ρR;
		w(1, i) = ρR*uR;
		w(2, i) = ρeR;
	}
	
	double tau;
	double λmax = 0;
	int counter = 0;
	for (double t = 0.0; t < T;)
	{
		for (int j = 0; j < NX; ++j)
		{
			if (u(w(0, j), w(1, j)) + c(w(0, j), w(2, j)) > λmax)
				λmax = u(w(0, j), w(1, j)) + c(w(0, j), w(2, j));
		}

		tau = std::min(CFL * h * 0.1 / λmax, 10e-2);
		t += tau;

		Eigen::Matrix3d Λ;
		Eigen::Matrix3d ΩT;
		Eigen::Matrix3d ΩT_Inv;
		Eigen::Matrix3d A;
		Eigen::Matrix3d ABS_Λ;

		for (int j = 1; j < NX-1; ++j)
		{
			double c2 = c(w(0, j), w(2, j)) * c(w(0, j), w(2, j));
			double u2 = u(w(0, j), w(1, j)) * u(w(0, j), w(1, j));
			double uj = u(w(0, j), w(1, j));
			double cj = c(w(0, j), w(2, j));

			Λ << uj + cj, 0.0, 0.0,
				0.0, uj, 0.0, 
				0.0, 0.0, uj - cj;

			ABS_Λ = Λ.cwiseAbs();

			ΩT << -uj * cj, cj, γ - 1.0,
				-c2, 0.0, γ - 1.0,
				uj* cj, -cj, γ - 1.0;
 
			ΩT_Inv << 1.0 / (2.0 * c2), -2.0 / (2.0 * c2), 1.0 / (2.0 * c2),
				(uj + cj) / (2.0 * c2), -2.0 * uj / (2.0 * c2), (uj - cj) / (2.0 * c2),
				1.0 / (2.0 * (γ - 1.0)), 0.0, 1.0 / (2.0 * (γ - 1.0));
			 
			A = ΩT_Inv * Λ * ΩT;

			w1.col(j) = w.col(j) - tau * A / (2.0 * h) * (w.col(j + 1) - w.col(j - 1)) +
				tau * (ΩT_Inv * ABS_Λ * ΩT) / (2.0 * h) * (w.col(j+1) - 2.0 * w.col(j) + w.col(j - 1));

		}
		w1.col(0) = w1.col(1);
		w1.col(NX-1) = w1.col(NX-2);
		w = w1;

	}
	std::ofstream Parameters;
	Parameters.open("dense.txt", std::ios_base::out | std::ios_base::trunc);
	for (int i = 0; i < NX; ++i)
	{
		Parameters << w(0, i) << std::endl;
	}
	Parameters.close();

	Parameters.open("pressure.txt", std::ios_base::out | std::ios_base::trunc);
	for (int i = 0; i < NX; ++i)
	{
		Parameters << P(w(2, i)) << std::endl;
	}
	Parameters.close();

	Parameters.open("dense2.txt", std::ios_base::out | std::ios_base::trunc);
	for (int i = 0; i < NX; ++i)
	{
		Parameters << u(w(0, i), w(1, i)) << std::endl;
	}
	Parameters.close();

	Parameters.open("energy.txt", std::ios_base::out | std::ios_base::trunc);
	for (int i = 0; i < NX; ++i)
	{
		Parameters << w(2, i)/w(0, i) << std::endl;
	}
	Parameters.close();
}


int main()
{
	Solver();

	return 0;
}
