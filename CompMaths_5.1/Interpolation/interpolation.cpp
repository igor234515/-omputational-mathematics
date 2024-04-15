#include <iostream>

double difference(double f1, double f2, double x1, double x2)
{
    return (f2-f1)/(x2-x1);
}

double Gorner_scheme(double x, int array_size, double *coordinates , double *functions)
{
    double ans = functions[array_size - 1];
    for(int i = array_size - 2; i >= 0; i--)
    {
        ans = ans*(x - coordinates[i]) + functions[i];
    }
    return ans;
}

int main()
{
    int len = 10;
    //std::cin >> len;
    double* x_coord = new double[len];
    double* results = new double[len];
    //double x_coord = { 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000 };
   
    for (int i = 0; i < len; ++i)
    {
        std::cin >> x_coord[i];

    }
    for (int i = 0; i < len; ++i)
    {
        std::cin >> results[i];
    }
    
    int counter = 1;

    for (int counter = 1; counter < len; ++counter)
    {
        for (int i = len - 1; i >= counter; i--)
        {
            results[i] = difference(results[i-1], results[i], x_coord[i-counter], x_coord[i]);
        }
    }

    double x_find = 2010;
    //std::cin >> x_find;
    std::cout << std::endl;
    /*for (int i = 0; i < len; ++i)
    {
        std::cout << results[i] << std::endl;
    }
    */
    std::cout << " Population = " << Gorner_scheme(x_find, len, x_coord, results);

    delete[] x_coord;
    delete[] results;
return 0;
}
/*
10
1910 1920 1930 1940 1950 1960 1970 1980 1990 2000
92228496
106021537
123202624
132164569
151325798
179323175
203211926
226545805
248709873
281421906
2010
*/
