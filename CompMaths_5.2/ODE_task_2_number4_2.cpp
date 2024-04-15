#include <iostream>
#include <cmath>

const double h = 0.001;
const double eps = 1e-6;

double func(double t, double y, double dy)
{
    return sqrt(1.0/t/t + 2.718281828*y*y/log(t) - exp(dy)*y);
}

void boundary_problem(double t_min, double t_max, double y_start, double y_end)
{
    double h_loc = h;
    int len = (t_max - t_min) / h_loc;
    double* t = new double[len];
    double* y = new double[len];
    double* dy = new double[len];
    double* ddy = new double[len];
    t[0] = t_min;
    y[0] = y_start;
    double delta1, delta2, delta;
    //сетка по времени
    for (int i = 1; i < len; ++i)
    {
        t[i] = t[i-1] + h_loc;
    }
    //начальнаое решение - прямая:
    double k = (y_end - y_start) / (t_max - t_min);
    dy[0] = k;
    ddy[0] = 0;
    for (int i = 1; i < len; ++i)
    {
        y[i] = y[0] + t[i]*k;
        dy[i] = k;
        ddy[0] = 0;
    }
    double k1, k2, k3, k4;
    /*
    for (int i = 1; i < len; ++i)
    {
        k1 = func(t[i-1], y[i-1], dy[i]);
        k2 = func(t[i-1] + h_loc/2.0, y[i-1] + h_loc/2.0 * k1, dy[i]);
        k3 = func(t[i-1] + h_loc/2.0, y[i-1] + h_loc/2.0 * k2, dy[i]);
        k4 = func(t[i], y[i-1] + h_loc * k3, dy[i]);

        y[i] = y[i-1] + h_loc / 6.0 *(k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }*/
    //std::cout << v[len-1].x << std::endl;
    /*
    delta1 = v[len-1].x - v_end;
    delta2 = delta1;
    double x1 = 0.0, x2;
    int counter = 1;
    while (delta1*delta2 > 0)
    {
        if (counter%4 == 0)
        {
            v[0].y = counter;
        }
        else
        {
            v[0].y = 1.0 / counter;
        }

        for (int i = 1; i < len; ++i)
        {
            k1 = func_sys(t[i-1], v[i-1]);
            k2 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k1);
            k3 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k2);
            k4 = func_sys(t[i], v[i-1] + h_loc * k3);

            v[i] = v[i-1] + h_loc / 6.0 *(k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        delta2 = v[len-1].x - v_end;

        counter *= 2;
    }
    x2 = v[0].y;

    if(x2 < x1)
    {
        x2 = x1;
        delta = delta2;
        delta2 = delta1;
        delta1 = delta;
        x1 = v[0].y;
    }

    while (fabs(delta1 / v_end) > eps)
    {
        v[0].y = (x2 + x1) / 2.0;

        for (int i = 1; i < len; ++i)
        {
            k1 = func_sys(t[i-1], v[i-1]);
            k2 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k1);
            k3 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k2);
            k4 = func_sys(t[i], v[i-1] + h_loc * k3);

            v[i] = v[i-1] + h_loc / 6.0 *(k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        delta = v[len-1].x - v_end;
        if (delta * delta1 > 0)
        {
            delta1 = delta;
            x1 = v[0].y;
        }
        else
        {
            delta2 = delta;
            x2 = v[0].y;
        }
        std::cout << x1 << " " << x2 << std::endl << delta1 << " " << delta2 << std::endl;
    }

    v[0].y = (x2 + x1) / 2.0;

    for (int i = 1; i < len; ++i)
    {
        k1 = func_sys(t[i-1], v[i-1]);
        k2 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k1);
        k3 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k2);
        k4 = func_sys(t[i], v[i-1] + h_loc * k3);

        v[i] = v[i-1] + h_loc / 6.0 *(k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    */
    for (int i = 0; i < len; ++i)
    {
        std::cout << y[i] << std::endl;
    }
    delete[] t;
    delete[] y;
    delete[] dy;
    delete[] ddy;
}

int main() {
    boundary_problem(2.718282, 7.389, 2.718282, 14.7781);
return 0;
}
