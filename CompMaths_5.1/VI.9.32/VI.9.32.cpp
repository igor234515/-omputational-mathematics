#include <iostream>
#include <cmath>

const double eps_step = 1e-7;
const double eps_func = 1e-7;

struct params
{
    double a;
};

double func(double coordinate, params p)
{
    return coordinate * exp(-1.0 * coordinate * coordinate) - p.a;
}

double difference(double coordinate, double step, params p)
{
    return (func(coordinate + step, p) - func(coordinate - step, p)) / (2.0 * step);
}

double solve(double x0, double step, params p)
{
    double x1;
    x1 = x0 - func(x0, p) / difference(x0, step, p);

    while (step > eps_step || fabs(x1 - x0) > eps_func)
    {
        x0 = x1;
        x1 = x0 - func(x0, p) / difference(x0, step, p);
        step = std::max(fabs(x1 - x0) / 4.0, eps_step / 10.0);
    }

    for (short i = 0; i < 5; ++i)
    {
        x0 = x1;
        x1 = x0 - func(x0, p) / difference(x0, step, p);
        step = std::max(fabs(x1 - x0) / 4.0, eps_step / 10.0);
    }
    return x1;
}

int main()
{
    double x0, step = 0.005;
    std::cout << "start point:" << std::endl;
    std::cin >> x0;
    params p;
    p.a = exp(-1.0 / 2.0) / sqrt(2.0) / 2.0;
    std::cout << "x = " << solve(x0, step, p) << std::endl;

    return 0;
}
