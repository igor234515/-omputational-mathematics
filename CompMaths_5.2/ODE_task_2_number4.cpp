#include <iostream>
#include <cmath>

const double h = 0.001;
const double eps = 1e-6;

struct vec_U{
    double x;
    double y;
};

vec_U func_sys(double t, vec_U u)
{
    vec_U ans;
    ans.x = u.y;
    ans.y = sqrt(1.0/t/t + 2.718281828*u.x*u.x/log(t) - exp(u.y)*u.x);
    return ans;
}

vec_U operator* (double t, vec_U u)
{
    vec_U ans;
    ans.x = t*u.x;
    ans.y = t*u.y;
    return ans;
}

vec_U operator+ (vec_U u1, vec_U u2)
{
    vec_U ans;
    ans.x = u1.x + u2.x;
    ans.y = u1.y + u2.y;
    return ans;
}

void boundary_problem(double t_min, double t_max, double v_start, double v_end)
{
    double h_loc = h;
    int len = (t_max - t_min) / h_loc;
    double* t = new double[len];
    vec_U* v = new vec_U[len];
    t[0] = t_min;
    v[0].x = v_start;
    v[0].y = 2.0067;
    double delta1, delta2, delta;

    for (int i = 1; i < len; ++i)
    {
        t[i] = t[i-1] + h_loc;
    }
    vec_U k1, k2, k3, k4;
    for (int i = 1; i < len; ++i)
    {
        k1 = func_sys(t[i-1], v[i-1]);
        k2 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k1);
        k3 = func_sys(t[i-1] + h_loc/2.0, v[i-1] + h_loc/2.0 * k2);
        k4 = func_sys(t[i], v[i-1] + h_loc * k3);

        v[i] = v[i-1] + h_loc / 6.0 *(k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
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
        std::cout << v[i].x << std::endl;
    }
    delete[] t;
    delete[] v;
}

int main() {
    boundary_problem(2.718282, 7.389, 2.718282, 14.7781);
return 0;
}
