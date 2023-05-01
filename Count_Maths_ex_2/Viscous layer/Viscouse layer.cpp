#include <iostream>
#include <fstream>

const double k = 1e-14;
const double φ = 0.2;
const double μ = 1e-3;
const double cf = 1e-4/101325.0;
const double p_sup0 = 12159000.0;
const double ρ0 = 1000.0;
//---------------------------
class Plast
{

public:
	int L = 500;
	int NX = 100;
	double T = 10 * 3600*24;
	double P0 = 10132500.0;
	double p_inj = 15198750.0;
	double p_prod = 5066250.0;
	double* p_layer = new double[NX];//for test, but it must be private
	double* d = new double[NX];
	double* c_diag = new double[NX];
	double* a_diag = new double[NX];
	double* b_diag = new double[NX];
	void print(double* ar);
	Plast();
	~Plast();
	void Solve(bool record, double time_output_init, double step_output);
private:

	//Grid
	double t = 0.5*3600;
	double h = L / NX;
	//For Thompson
	double* p = new double[NX];
	double* q = new double[NX];

	//Function prototypes
	double ρ(double p);
	double ρ_minus_1_2(int i);
	double ρ_plus_1_2(int i);
	void makematrix();
	void Thompson();
};
//--------------------------------
Plast::Plast()
{
	p_layer[0] = p_inj;
	p_layer[NX - 1] = p_prod;
	for (uint8_t i = 1; i < NX - 1; i++)
	{
		p_layer[i] = P0;
	}
	std::cout << "Init Values" << std::endl;;
	print(p_layer);
}

Plast::~Plast()
{
	delete[] p_layer;
	delete[] a_diag;
	delete[] b_diag;
	delete[] c_diag;
	delete[] d;
	delete[] p;
	delete[] q;
	std::cout << "Object deleted";
}

double Plast::ρ(double p)
{
	return ρ0 * (1.0 + cf * (p - p_sup0));
}

double Plast::ρ_plus_1_2(int i)
{
	double ro_1_2;
	if (p_layer[i] >= p_layer[i + 1]) ro_1_2 = ρ(p_layer[i]);
	else ro_1_2 = ρ(p_layer[i + 1]);
	return ro_1_2;
}

double Plast::ρ_minus_1_2(int i)
{
	double ro_1_2;
	if (p_layer[i - 1] >= p_layer[i]) ro_1_2 = ρ(p_layer[i - 1]);
	else ro_1_2 = ρ(p_layer[i]);
	return ro_1_2;
}

void Plast::makematrix()
{
	d[0] = p_inj;
	d[NX - 1] = p_prod;
	a_diag[0] = 1;
	a_diag[NX - 1] = 1;
	b_diag[0] = 0;
	b_diag[NX - 1] = 0;
	c_diag[0] = 0;
	c_diag[NX - 1] = 0;

	for (uint8_t i = 1; i < NX - 1; ++i)
	{
		c_diag[i] = k * ρ_minus_1_2(i) / (μ * h * h);
		b_diag[i] = k * ρ_plus_1_2(i) / (μ * h * h);
		a_diag[i] = -c_diag[i] - b_diag[i] - φ * cf * ρ0 / t;
		d[i] = -φ * cf * ρ0 / t * p_layer[i];
	}
}

void Plast::Thompson()
{
	//p[1] = -b_diag[0] / a_diag[0]; //In general
	//q[1] = d[0] / a_diag[0];
	p[1] = 0;
	q[1] = d[0];
	for (uint8_t i = 1; i < NX-1; ++i)
	{
		p[i + 1] = -b_diag[i] / (c_diag[i] * p[i] + a_diag[i]);
		q[i + 1] = (d[i] - c_diag[i] * q[i]) / (c_diag[i] * p[i] + a_diag[i]);
	}
	p_layer[NX - 1] = (d[NX - 1] - c_diag[NX - 1] * q[NX - 1]) / (c_diag[NX - 1] * p[NX - 1] + a_diag[NX - 1]);
	for (int i = NX - 2; i >= 0  ; --i)
	{
		p_layer[i] = p_layer[i + 1] * p[i + 1] + q[i + 1];
	}
}

void Plast::Solve(bool record = false, double time_output_init = 0,  double step_output = 0.1)
{
	double time = 0;
	double time_output = time_output_init*3600*24;//In days
	double dt_out_results = step_output*3600*24;
	if (record == false)
	{
		while (time <= T)
		{
			makematrix();
			Thompson();
			time += t;
			if (time >= time_output)
			{
				std::cout << "Pressure distribution" << std::endl;
				print(p_layer);
				std::cout << std::endl;
				time_output += dt_out_results;
			}

		}
	}
	else
	{
		std::ofstream Pressure;
		Pressure.open("Pressure.txt", std::ios_base::out | std::ios_base::trunc);
		while (time <= T)
		{
			makematrix();
			Thompson();
			time += t;
			if (time >= time_output)
			{
				for (uint8_t i = 0; i < NX; ++i)
				{
					Pressure << p_layer[i] << std::endl;
				}
				time_output += dt_out_results;
			}

		}
		Pressure.close();
	}
}

void Plast::print(double* ar)
{
	for (int i = 0; i < NX; ++i)
	{
		std::cout << ar[i] << " ";
	}
}

//--------------------

int main()
{
	Plast oil;
	oil.Solve(true);
	return 0;
}