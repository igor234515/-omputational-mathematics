#include <iostream>
#include <cmath>
#include <fstream>

const long double mu = 0.012277471, nu = 1.0 - 0.012277471;

struct params{
    long double x;
    long double y;
    long double u;
    long double v;
};

params operator* (long double t, params a)
{
    params ans;
    ans.x = t*a.x;
    ans.y = t*a.y;
    ans.u = t*a.u;
    ans.v = t*a.v;
    return ans;
}

params operator+ (params u1, params u2)
{
    params ans;
    ans.x = u1.x + u2.x;
    ans.y = u1.y + u2.y;
    ans.u = u1.u + u2.u;
    ans.v = u1.v + u2.v;
    return ans;
}

params system_diff(long double t, params s)
{
    params ans;
    long double A, B;
    A = sqrt((s.x + mu)*(s.x + mu) + s.y*s.y)*((s.x + mu)*(s.x + mu) + s.y*s.y);
    B = sqrt((s.x - nu)*(s.x - nu) + s.y*s.y)*((s.x - nu)*(s.x - nu) + s.y*s.y);
    ans.x = s.u;
    ans.u = s.x + 2.0*s.v - nu*(s.x + mu)/A - mu*(s.x - nu)/B;
    ans.y = s.v;
    ans.v = s.y - 2.0*s.u - nu*s.y/A - mu*s.y/B;
    return ans;
}

int main()
{
    int n = 7;
    const long double tau = 1e-4;
    const long double T = 17.0652165601579625588917206249;
    long double t_start = 0.0;
    long double t_end = 100.0*T;
    const long long int len = (T - t_start) / tau;
    long double b[7] = {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0};
    long double c[7] = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
    //для одного периода
    params* state = new params[len];
    state[0].x = 0.994;
    state[0].y = 0.0;
    state[0].u = 0.0;
    state[0].v = -2.00158510637908252240537862224;

    const long double butcher_matrix[7][7] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {3.0/40.0, 9.0/40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {44.0/45.0, -56.0/15.0, 32.0/9.0, 0.0, 0.0, 0.0, 0.0},
    {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0.0, 0.0, 0.0},
    {9017.0/3168.0,	-355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0, 0.0, 0.0},
    {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0} };

    for (unsigned int i = 1; i < len; ++i)
    {
        params k[7];
        params mult;
        k[0] = system_diff(tau*i, state[i-1]);
        for (short j = 1; j < 7; ++j)
        {
            mult.x = 0.0;
            mult.y = 0.0;
            mult.u = 0.0;
            mult.v = 0.0;
            for (short p = 0; p < j-1; ++p)
            {
                mult = mult + butcher_matrix[j][p]*k[p];
            }
            k[j] = system_diff(tau*i + tau*j, state[i-1] + tau*mult);
        }
        mult.x = 0.0;
        mult.y = 0.0;
        mult.u = 0.0;
        mult.v = 0.0;
        for (short j = 0; j < 7; ++j)
        {
            mult = mult + b[j]*k[j];
        }
        state[i] = state[i-1] + tau*mult;
    }

   // for (unsigned int i = 0; i < len/10000; ++i)
   // {
    //    std::cout << "step " << i*10000 << " : " << state[i*10000].x << " " << state[i*10000].y << std::endl;
   // }

    std::ofstream outfileX;
    std::ofstream outfileY;
    outfileX.open("Trajectory_X.txt", std::ios_base::out | std::ios_base::trunc);
    outfileY.open("Trajectory_Y.txt", std::ios_base::out | std::ios_base::trunc); //std::ios_base::app
    for (int i = 0; i < len; ++i)
    {
        
        outfileX << state[i].x << std::endl;
        outfileY << state[i].y << std::endl;
    }
    outfileX.close();
    outfileY.close();

    delete[] state;

return 0;
}
