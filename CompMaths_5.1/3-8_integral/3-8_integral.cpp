#include <iostream>

int main()
{
    int len;
    double step;
    std::cout << "count of point; step width" << std::endl;
    std::cin >> len >> step;
    double* results = new double[len];

    for (int i = 0; i < len; ++i)
    {
        std::cin >> results[i];
    }

    double integral = 0;

    if ((len-1)%3 == 0)
    {
        for (int i = 1; i < len-1; ++i)
        {
            if (i%3 == 0)
            {
                integral += 2.0 * results[i];
            }
            else
            {
                integral += 3.0 * results[i];
            }
        }
        integral += results[0] + results[len-1];
        integral *= 3.0 * step / 8.0;
    }

    else if ((len-1)%3 == 1)
    {
        for (int i = 1; i < len-2; ++i)
        {
            if (i%3 == 0)
            {
                integral += 2.0 * results[i];
            }
            else
            {
                integral += 3.0 * results[i];
            }
        }
        integral += results[0] + results[len-2];
        integral *= 3.0 * step / 8.0;
        integral += step * (results[len-1] + results[len-2]) / 2.0;
    }

    else
    {
        for (int i = 1; i < len-3; ++i)
        {
            if (i%3 == 0)
            {
                integral += 2.0 * results[i];
            }
            else
            {
                integral += 3.0 * results[i];
            }
        }
        integral += results[0] + results[len-3];
        integral *= 3.0 * step / 8.0;
        integral += step * (results[len-1] + 4.0*results[len-2] + results[len-3]) / 3.0;
    }

    std::cout << integral;
return 0;
}
