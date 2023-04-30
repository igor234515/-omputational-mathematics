#include <iostream>
#include <Eigen/Dense>



//Объект - пласт с заданными параметрами

class Plast
{
public:
	//Public variables
	int L = 500;
	int NX = 100;
	double T = 0.1 * 24;
	double P0 = 10132500;
	double p_inj = 15198750;
	double p_prod = 5066250;

private:
	//Private variables
	double k = 1e-14;
	double φ = 0.2;
	double μ = 1e-3;
	double cf = 1e-4;
	double p_sup0 = 12159000;
	double ρ_sup = 1000;
	//Grid
	
	double t = 0.5;
	int Nt = T / t;
	double h = L / NX;
	//Function prototypes
public:

	//Constructor - creation a layer with start conditions - or not?
	Plast()
	{
		p_layer[0] = p_inj;
		p_layer[NX - 1] = p_prod;
			for (int i = 1; i < NX-1; i++)
			{
				p_layer[i] = P0;
	}
	//Function prototypes

};


int main()
{
	Plast oil;
	for (int i = 0; i < oil.NX; i++)
	{
		std::cout << oil.p_layer[i] << "\t";
	}
	return 0;
}