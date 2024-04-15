#include <iostream>
#include <cmath>

const double h = 0.001;

struct vec_U{
    double x;
    double y;
};

double func(double x, double y)
{
    return y;
}

vec_U func_sys(double t, vec_U u)
{
    vec_U ans;
    ans.x = u.y;
    ans.y = cos(3.0*t) - 4.0 * u.x;
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

void forward_Euler_method()
{
    double x_min, x_max;
    std::cout << "borders of x coordinate" << std::endl;
    std::cin >> x_min >> x_max;
    int len = (x_max - x_min) / h;
    double* x_coord = new double[len];
    double* y = new double[len];
    x_coord[0] = x_min;
    std::cout << "boundary conditions: y (x start)" << std::endl;
    std::cin >> y[0];
    for (int i = 1; i < len; ++i)
    {
        x_coord[i] = x_coord[i-1] + h;
    }
    for (int i = 1; i < len; ++i)
    {
        y[i] = y[i-1] + h*func(x_coord[i-1], y[i-1]);
    }

    for (int i = 0; i < len; ++i)
    {
        std::cout << y[i] << std::endl;
    }
    delete[] x_coord;
    delete[] y;
}

void not_forward_Euler_method()
{
    double x_min, x_max;
    std::cout << "borders of x coordinate" << std::endl;
    std::cin >> x_min >> x_max;
    int len = (x_max - x_min) / h;
    double* x_coord = new double[len];
    double* y = new double[len];
    x_coord[0] = x_min;
    std::cout << "boundary conditions: y (x start)" << std::endl;
    std::cin >> y[0];
    for (int i = 1; i < len; ++i)
    {
        x_coord[i] = x_coord[i-1] + h;
    }
    for (int i = 1; i < len; ++i)
    {
        y[i] = y[i-1] + h*func(x_coord[i], y[i-1]);
    }

    for (int i = 0; i < len; ++i)
    {
        std::cout << y[i] << std::endl;
    }
    delete[] x_coord;
    delete[] y;
}

void center_scheme_method()
{
    double x_min, x_max;
    std::cout << "borders of x coordinate" << std::endl;
    std::cin >> x_min >> x_max;
    int len = (x_max - x_min) / h;
    double* x_coord = new double[len];
    double* y = new double[len];
    x_coord[0] = x_min;
    x_coord[1] = x_min + h;
    std::cout << "boundary conditions: y (x start), y (x start + h)" << std::endl;
    std::cin >> y[0]>> y[1];
    for (int i = 2; i < len; ++i)
    {
        x_coord[i] = x_coord[i-1] + h;
    }
    for (int i = 2; i < len; ++i)
    {
        y[i] = y[i-2] + 2.0*h*func(x_coord[i-1], y[i-1]);
    }

    for (int i = 0; i < len; ++i)
    {
        std::cout << y[i] << std::endl;
    }
    delete[] x_coord;
    delete[] y;
}

void Euler_method_with_recalculation()
{
    double x_min, x_max;
    std::cout << "borders of x coordinate" << std::endl;
    std::cin >> x_min >> x_max;
    int len = (x_max - x_min) / h;
    double* x_coord = new double[len];
    double* y = new double[len];
    x_coord[0] = x_min;
    std::cout << "boundary conditions: y (x start)" << std::endl;
    std::cin >> y[0];
    for (int i = 1; i < len; ++i)
    {
        x_coord[i] = x_coord[i-1] + h;
    }
    double k;
    for (int i = 1; i < len; ++i)
    {
        k = func(x_coord[i-1], y[i-1]);
        y[i] = y[i-1] + h*func(x_coord[i-1] + h/2.0, y[i-1] + k*h/2.0);
    }

    for (int i = 0; i < len; ++i)
    {
        std::cout << y[i] << std::endl;
    }
    delete[] x_coord;
    delete[] y;
}

void R_K_method_4_order()
{
    double x_min, x_max;
    std::cout << "borders of x coordinate" << std::endl;
    std::cin >> x_min >> x_max;
    int len = (x_max - x_min) / h;
    double* x_coord = new double[len];
    double* y = new double[len];
    x_coord[0] = x_min;
    std::cout << "boundary conditions: y (x start)" << std::endl;
    std::cin >> y[0];
    for (int i = 1; i < len; ++i)
    {
        x_coord[i] = x_coord[i-1] + h;
    }
    double k1, k2, k3, k4;
    for (int i = 1; i < len; ++i)
    {
        k1 = func(x_coord[i-1], y[i-1]);
        k2 = func(x_coord[i-1] + h/2.0, y[i-1] + k1*h/2.0);
        k3 = func(x_coord[i-1] + h/2.0, y[i-1] + k2*h/2.0);
        k4 = func(x_coord[i], y[i-1] + k3*h);

        y[i] = y[i-1] + h*(k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0);
    }

    for (int i = 0; i < len; ++i)
    {
        std::cout << y[i] << std::endl;
    }
    delete[] x_coord;
    delete[] y;
}

void R_K_method_4_order_system(double h_loc)
{
    double t_min, t_max;
    std::cout << "borders of t coordinate" << std::endl;
    std::cin >> t_min >> t_max;
    int len = (t_max - t_min) / h_loc;
    double* t = new double[len];
    vec_U* v = new vec_U[len];
    t[0] = t_min;
    std::cout << "boundary conditions U start:" << std::endl;
    std::cin >> v[0].x >> v[0].y;
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

    for (int i = 0; i < len; ++i)
    {
        std::cout << v[i].x << " " << v[i].y << std::endl;
    }
    delete[] t;
    delete[] v;
}

double g_x(double x)
{
    return x*x - 3.0;
}

double h_x(double x)
{
    return (x*x - 3.0) * cos(x);
}

double f_x(double x)
{
    return 2.0 - 6.0*x + 2.0*x*x*x + (x*x - 3.0)*exp(x)*sin(x)*(1.0 + cos(x)) + cos(x)*(exp(x) + x*x - 1.0 +  x*x*x*x - 3.0*x*x);
}

//y'' + g(x)*y' + h(x)*y = f(x)
void boundary_problem(double h_loc)
{
    const double eps = 1e-4;
    int len;
    double t_start, t_end, y_last;
    std::cout << "bounds: start, end" << std::endl;
    std::cin >> t_start >> t_end;
    len = (t_end - t_start) / h_loc + 1;
    double* t = new double[len];
    t[0] = t_start;
    for (int i = 1; i < len; ++i)
    {
        t[i] = t[i-1] + h_loc;
    }
    double* a = new double[len];
    double* b = new double[len];
    double* c = new double[len];
    double* d = new double[len];
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = 0.0;
    a[len-1] = 0.0;
    b[len-1] = 1.0;
    c[len-1] = 0.0;
    for (int i = 1; i < len-1; ++i)
    {
        a[i] = 1.0 + g_x(t[i])*h_loc;
        b[i] = -2.0 - g_x(t[i])*h_loc + h_x(t[i])*h_loc*h_loc;
        c[i] = 1.0;
        d[i] = f_x(t[i])*h_loc*h_loc;
    }

    double* y = new double[len];
    std::cout << "boundary conditionals:" << std::endl;
    std::cin >> y[0] >> y_last;

    //подбор y[1]. 2 перебираемзначения, пока не найдём такие, при которых отклонения от граничного условия имеют разные знаки.
    double delta1, delta2, delta;
    y[1] = y[0];

    for (int i = 2; i < len; ++i)
    {
        y[i] = (h_loc*h_loc*f_x(t[i-1]) - y[i-2] - y[i-1]*b[i-1]) / a[i-1];
    }
    delta1 = y[len-1] - y_last;
    int counter = 1;
    delta2 = delta1;
    while (delta1*delta2 > 0)
    {
        counter *= 2;
        if (y[0] == 0.0)
        {
            y[1] = counter * (-1.0)*((counter % 4)/2 + 1);
        }
        else
        {
            y[1] = y[0] * counter * (-1.0)*((counter % 4)/2 + 1);
        }

        for (int i = 2; i < len; ++i)
        {
            y[i] = (h_loc*h_loc*f_x(t[i-1]) - y[i-2] - y[i-1]*b[i-1]) / a[i-1];
        }
        delta2 = y[len-1] - y_last;
    }
    double y1, y2, y_new;
    if ((y[0] - y[1]) > 0)
    {
        y2 = y[0];
        y1 = y[1];
        delta = delta1;
        delta1 = delta2;
        delta2 = delta;
    }
    else
    {
        y2 = y[1];
        y1 = y[0];
    }
    while (fabs(delta1 / y_last) > eps)
    {
        y_new = (y2 + y1) / 2.0;
        y[1] = y_new;
        for (int i = 2; i < len; ++i)
        {
            y[i] = (h_loc*h_loc*f_x(t[i-1]) - y[i-2] - y[i-1]*b[i-1]) / a[i-1];
        }
        delta = y[len-1] - y_last;
        if (delta*delta2 > 0)
        {
            y2 = y_new;
            delta2 = delta;
        }
        else
        {
            y1 = y_new;
            delta1 = delta;
        }
    }
    y[1] = (y1 + y2) / 2.0;
    for (int i = 2; i < len; ++i)
    {
        y[i] = (h_loc*h_loc*f_x(t[i-1]) - y[i-2] - y[i-1]*b[i-1]) / a[i-1];
    }

    for (int i = 0; i < len; ++i)
    {
        std::cout << y[i] << std::endl;
    }

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] t;
    delete[] y;
}

int main() {
    boundary_problem(h);
    //R_K_method_4_order_system(h);
return 0;
}
