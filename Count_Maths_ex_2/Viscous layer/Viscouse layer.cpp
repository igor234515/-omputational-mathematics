#include <iostream>

const double k = 1e-14;
const double φ = 0.2;
const double μ = 1e-3;
const double cf = 1e-4;
const double p_sup0 = 12159000;
const double ρ0 = 1000;
//---------------------------
class Plast
{

public:
	int L = 500;
	int NX = 100;
	double T = 0.1 * 24;
	double P0 = 10132500;
	double p_inj = 15198750;
	double p_prod = 5066250;
	double* p_layer = new double[NX];//for test, but it must be private
	double* d = new double[NX];
	double* c_diag = new double[NX];
	double* a_diag = new double[NX];
	double* b_diag = new double[NX];
	void print(double* ar);
	Plast();
	~Plast();
private:

	//Grid
	double t = 0.5;
	int Nt = T / t;
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
	makematrix();
}

double Plast::ρ(double p)
{
	return ρ0 * (1 + cf * (p - p_sup0));
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
	std::cout << "a of matrix" << std::endl;
	oil.print(oil.a_diag);
	return 0;
}