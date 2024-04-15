#include <iostream>
#include <cmath>

double func (double x)
{
    return sin(100.0*x) * exp(-1.0 * x*x) * cos(2.0 * x);
}

double integrate (int len, double a, double b)
{
    double step = (b - a) / len;
    double* results = new double[len];

    for (int i = 0; i < len; ++i)
    {
        results[i] = func(a + i * step);
    }

    double integral = 0;

    for (int i = 1; i < len-1 - (len-1)%2; ++i)
    {
        integral += 2.0 * (i%2 + 1.0) * results[i];
    }
    integral += results[0] + results[len - 1 - (len-1)%2];
    integral *= step / 3.0;

    if (len%2 == 0)
    {
        integral += step*(results[len-1] + results[len-2])/2;
    }
    delete[] results;

    return integral;
}

int main()
{
    //int len;
    double did = 15000;
    double a = 0.0, b = 0.8;
    /*std::cout << "count of point; step width" << std::endl;
    std::cin >> len >> step;*/
    //std::cout << "count of point; start point; end point" << std::endl;
    //std::cin >> len >> a >> b;
    /*for (int i = 1000; i < 1500000; i *= 2.0)
    {
        std::cout << "count of point = " << i << " : " << integrate(i, a, b) << std::endl;
    }
    */
    std::cout << " I = " << integrate(did, a, b) << std::endl;
return 0;
}
