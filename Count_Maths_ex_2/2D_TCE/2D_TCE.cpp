#include <iostream>
#include <cmath>
const double PI = 3.14;

class TC2D
{
public:
	TC2D();
	~TC2D();
	void print_plot();
	void print_layer();
	double λ = 1e-4;
	int X = 1;
	double Y = 1;

private:
	int NX = 10;
	int NY = 10;
	double h_x = X / NX;
	double h_y = Y / NY;
	double** xy_plot = new double* [NX];
};

TC2D::TC2D()
{
	for (int i = 0; i < NX; ++i)
	{
		xy_plot[i] = new double[NY];
	}
	for (int i = 0; i < NX; ++i)
	{
		for (int j = 0; j < NY-1; ++j)
		{
			xy_plot[i][j] = cos(PI * h_x * i) * sin(5 * PI * h_y * j);
		}
		xy_plot[i][NY - 1] = xy_plot[i][0];
	}
}

TC2D::~TC2D()
{
	delete[] xy_plot;
	std::cout << "Object deleted";
}

void TC2D::print_plot()
{
	for (int i = 0; i < NX; ++i)
	{
		for (int j = 0; j < NY; ++j)
		{
			std::cout << xy_plot[i][j] << " ";
		}
	}
}
void TC2D::print_layer()
{
	for (int i = 0; i < NY; ++i)
	{
		std::cout << xy_plot[0][i]<<" ";
	}
}
int main()
{
	TC2D eq;
	eq.print_layer();
	return 0;
}