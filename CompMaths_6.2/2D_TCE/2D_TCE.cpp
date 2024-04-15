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
	void makematrix(int j);
	void Thompson(double* a_diag, double* b_diag, double* c_diag, double* d, int N, double* return_aray);
	void Plot_Solve();
	double λ = 1e-4;
	double a = 25 * λ;
	int X = 1;
	double Y = 1;
	int T = 1;
	const double CFL = 0.1;
	double** xy_plot = new double* [NX];

private:
	int NX = 10;
	int NY = 10;
	double h_x = X / NX;
	double h_y = Y / NY;
	double t = CFL*h_x*h_x/a;
	int NT = T / t;
	double time = 0;
	
	//For thompson
	double* a_diag_x = new double[NX];
	double* b_diag_x = new double[NX];
	double* c_diag_x = new double[NX];
	double* d_x = new double[NX];
	double* a_diag_y = new double[NY];
	double* b_diag_y = new double[NY];
	double* c_diag_y = new double[NY];
	double* d_y = new double[NY];
};

TC2D::TC2D()
{
	for (int i = 0; i < NX; ++i)
	{
		xy_plot[i] = new double[NY];
	}
	for (int j = 0; j < NY; ++j)
	{
		xy_plot[0][j] = sin(5.0 * PI * h_y * j) * exp(-50.0 * PI * PI * λ * time);
		xy_plot[NX - 1][j] = -sin(5.0 * PI * h_y * j) * exp(-50.0 * PI * PI * λ * time);
	}

	for (int i = 1; i < NX-1; ++i)
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
	for (int i = 0; i < NX; ++i)
	{
		delete[] xy_plot[i];
	}
	delete[] a_diag_x;
	delete[] b_diag_x;
	delete[] c_diag_x;
	delete[] d_x;
	delete[] a_diag_y;
	delete[] b_diag_y;
	delete[] c_diag_y;
	delete[] d_y;
	std::cout << "Object deleted";
}
void TC2D::makematrix(int j)
{
	a_diag_x[0] = 1;
	b_diag_x[0] = 0;
	c_diag_x[0] = 0;
	d_x[0] = 0;
	d_x[NX - 1] = 0;
	for (int i = 1; i < NX-1; ++i)
	{
		a_diag_x[i] = 1 - 2* CFL;
		b_diag_x[i] = -CFL;
		c_diag_x[i] = -CFL;
		d_x[i] = xy_plot[i][j];
	}
	c_diag_x[NX - 1] = 0;
	a_diag_x[NX - 1] = 1;
	b_diag_x[NX - 1] = 0;
}
void TC2D::Thompson(double* a_diag, double* b_diag, double* c_diag, double* d, int N, double* return_aray)
{
	double* p = new double[N];
	double* q = new double[N];
	p[1] = -b_diag[0] / a_diag[0];
	q[1] = d[0] / a_diag[0];
	for (uint8_t i = 1; i < N - 1; ++i)
	{
		p[i + 1] = -b_diag[i] / (c_diag[i] * p[i] + a_diag[i]);
		q[i + 1] = (d[i] - c_diag[i] * q[i]) / (c_diag[i] * p[i] + a_diag[i]);
	}
	return_aray[N - 1] = (d[N - 1] - c_diag[N - 1] * q[N - 1]) / (c_diag[N - 1] * p[N - 1] + a_diag[N - 1]);
	for (int i = N - 2; i >= 0; --i)
	{
		return_aray[i] = return_aray[i + 1] * p[i + 1] + q[i + 1];
	}
}
void TC2D::Plot_Solve()
{
	for (int j = 0; j < NY; ++j)
	{
		makematrix(j);

	}
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
	std:: cout << std::endl;
}

int main()
{
	TC2D equation;
	equation.print_layer();
	return 0;
}