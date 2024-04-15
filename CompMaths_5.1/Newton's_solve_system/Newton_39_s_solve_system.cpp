#include <iostream>
#include <cmath>

const double eps_step = 1e-7;
const double eps_func = 1e-7;

struct point
{
    double x;
    double y;
};

struct matrix_2
{
    // number of line _ number of column
    double first_first;
    double first_second;
    double second_first;
    double second_second;
};

point func(point u)
{
    point ans;
    ans.x = u.x * u.x + u.y * u.y - 1.0;
    ans.y = tan(u.x)-u.y;
    return ans;
}

double det_2 (matrix_2 A)
{
    return A.first_first * A.second_second - A.first_second * A.second_first;
}

matrix_2 J(point u, double step)
{
    matrix_2 J_matrix;

    point u_plus = u, u_minus = u;
    u_plus.x += step;
    u_minus.x -= step;
    J_matrix.first_first = (func(u_plus).x - func(u_minus).x) / (2.0 * step);
    J_matrix.second_first = (func(u_plus).y - func(u_minus).y) / (2.0 * step);
    u_plus.x = u.x;
    u_minus.x = u.x;
    u_plus.y += step;
    u_minus.y -= step;
    J_matrix.first_second = (func(u_plus).x - func(u_minus).x) / (2.0 * step);
    J_matrix.second_second = (func(u_plus).y - func(u_minus).y) / (2.0 * step);

    return J_matrix;
}

matrix_2 reverse_matrix (matrix_2 A)
{
    double det = det_2(A);
    matrix_2 help;
    help.first_first = A.second_second / det;
    help.second_second = A.first_first / det;
    help.first_second = (-1.0) * A.first_second / det;
    help.second_first = (-1.0) * A.second_first / det;
    return help;
}

point multiplication (matrix_2 A, point v)
{
    point ans;
    ans.x = A.first_first * v.x + A.first_second * v.y;
    ans.y = A.second_first * v.x + A.second_second * v.y;
    return ans;
}

point solve(point u0, double step0)
{
    double step;
    point u1;
    matrix_2 Jacobian = J(u0, step0);
    Jacobian = reverse_matrix(Jacobian);
    point delta = multiplication(Jacobian, func(u0));
    u1.x = u0.x - delta.x;
    u1.y = u0.y - delta.y;
    step = std::max( std::max(fabs(u1.x - u0.x)/4.0, fabs(u1.y - u0.y)/4.0), eps_step/10.0);

    while (step > eps_step || std::max(fabs(u1.x - u0.x), fabs(u1.y - u0.y)) > eps_func)
    {
        u0 = u1;
        Jacobian = J(u0, step);
        Jacobian = reverse_matrix(Jacobian);
        point delta = multiplication(Jacobian, func(u0));
        u1.x = u0.x - delta.x;
        u1.y = u0.y - delta.y;
        step = std::max( std::max(fabs(u1.x - u0.x)/4.0, fabs(u1.y - u0.y)/4.0), eps_step/10.0);
    }

    for (short i = 0; i < 5; ++i)
    {
        u0 = u1;
        Jacobian = J(u0, step);
        Jacobian = reverse_matrix(Jacobian);
        point delta = multiplication(Jacobian, func(u0));
        u1.x = u0.x - delta.x;
        u1.y = u0.y - delta.y;
        step = std::max( std::max(fabs(u1.x - u0.x)/4.0, fabs(u1.y - u0.y)/4.0), eps_step/10.0);
    }
    return u1;
}

int main()
{
    point u0;
    double step0;
    //std::cout << "start point:" << std::endl;
    //std::cin >> u0.x >> u0.y;
    //std::cout << "start step:" << std::endl;
    //std::cin >> step0;
    u0.x = 0.65;
    u0.y = 0.76;
    step0 = 0.4;
    point ans = solve(u0, step0);
    std::cout << "x1 = " << ans.x << " " << "y1 = " << ans.y << std::endl;

    u0.x = -0.9;
    u0.y = -0.4;
    step0 = 0.05;
    ans = solve(u0, step0);
    std::cout << "x2 = " << ans.x << " " << "y2 = " << ans.y << std::endl;
return 0;
}
