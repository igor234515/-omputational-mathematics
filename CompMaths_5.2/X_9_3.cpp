#include <iostream>
#include <cmath>
#include <fstream>

const double h = 0.001;
const double eps = 1e-4;

struct vec_U{
    double x;
    double y;
};

vec_U func_sys(double t, vec_U u)
{
    vec_U ans;
    ans.x = u.y;
    ans.y = 1000.0*(1.0 - u.y*u.y)*u.y - u.x;
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

vec_U boundary_problem(double t_min, double t_max, double vx_start, double vy_start)
{
    std::ofstream outfile;

    double h_loc = h;
    int len = (t_max - t_min) / h_loc + 1;
    double* t = new double[len];
    vec_U* v = new vec_U[len];
    t[0] = t_min;
    v[0].x = vx_start;
    v[0].y = vy_start;

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

    outfile.open("X.9.3.txt - new", std::ios_base::out | std::ios_base::trunc);
    for (int i = 0; i < len; ++i)
    {
        outfile << v[i].x << std::endl;
    }
    outfile.close();

    vec_U ans;
    ans.x = v[len-1].x;
    ans.y = v[len-1].y;
    delete[] t;
    delete[] v;
    return ans;
}

int main() {
    vec_U param;
    param.x = 0.0;
    param.y = 0.001;
    double T0 = 0.0, T = 1000.0;

 /*   for (double t = T0; t < T; t += 1.0)
    {
        std::cout << "modelling from " << t << " to " << t + 1.0 << std::endl;
        param = boundary_problem(t, t + 1.0, param.x, param.y);
    } */
return 0;
}
