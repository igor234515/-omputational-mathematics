#include <iostream>
#include <cmath>
#include <fstream>

double L = 10.0;
double gamma = 5.0 / 3.0;
double T = 0.015;
double tau = 0.0001;

double v_l = 0.0;
double rho_l = 13.0;
double p_l = 1010000.0;

double v_r = 0.0;
double rho_r = 1.3;
double p_r = 101000.0;

double CFL = 0.01;
int NX = 100;
double h = 2.0 * L / NX;

struct state {
    double rho;
    double rho_u;
    double rho_e;
};

double pressure (state U) {
    return (gamma-1.0)*U.rho_e;
}

double find_c (state U) {
    return sqrt(gamma * pressure(U) / U.rho);
}

void omega_U (state U, double **adress) {
    double present_c = find_c(U);
    adress[0][0] = -1.0 * U.rho_u * present_c / U.rho;
    adress[0][1] = present_c;
    adress[1][0] = -1.0 * present_c * present_c;
    adress[1][1] = 0.0;
    adress[2][0] = U.rho_u * present_c / U.rho;
    adress[2][1] = -1.0 * present_c;
    for (int i = 0; i < 3; ++i)
    {
        adress[i][2] = gamma - 1.0;
    }
}

void reverse_omega_U (state U, double **adress) {
    double present_c = find_c(U);
    adress[0][0] = 0.5 / present_c / present_c;
    adress[0][1] = -1.0 / present_c / present_c;
    adress[0][2] = 0.5 / present_c / present_c;
    adress[1][0] = 0.5 * (U.rho_u / U.rho + present_c) / present_c / present_c;
    adress[1][1] = -1.0 * (U.rho_u / U.rho) / present_c / present_c;
    adress[1][2] = 0.5 * (U.rho_u / U.rho - present_c) / present_c / present_c;
    adress[2][0] = 0.5 / (gamma - 1.0);
    adress[2][1] = 0.0;
    adress[2][2] = 0.5 / (gamma - 1.0);
}

void lambda_U (state U, double **adress) {
    double present_c = find_c(U);
    adress[0][0] = U.rho_u / U.rho + present_c;
    adress[0][1] = 0.0;
    adress[0][2] = 0.0;
    adress[1][0] = 0.0;
    adress[1][1] = U.rho_u / U.rho;
    adress[1][2] = 0.0;
    adress[2][0] = 0.0;
    adress[2][1] = 0.0;
    adress[2][2] = U.rho_u / U.rho - present_c;
}

void mod_lambda_U (state U, double **adress) {
    double present_c = find_c(U);
    adress[0][0] = fabs(U.rho_u / U.rho + present_c);
    adress[0][1] = 0.0;
    adress[0][2] = 0.0;
    adress[1][0] = 0.0;
    adress[1][1] = fabs(U.rho_u / U.rho);
    adress[1][2] = 0.0;
    adress[2][0] = 0.0;
    adress[2][1] = 0.0;
    adress[2][2] = fabs(U.rho_u / U.rho - present_c);
}

double mult_lines(int length, double *str, double *column) {
    double ans = 0.0;
    for (int i = 0; i < length; ++i)
    {
        ans += str[i] * column[i];
    }
    return ans;
}

void transponse(int len, double **matrix) {
    double var;
    for (int i = 0; i < len - 1; ++i)
    {
        for (int j = i + 1; j < len; ++j)
        {
            var = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = var;
        }
    }
}

void mult_matrix(int len, double **first, double **second) {
    double** array1 = new double*[len];
    for (int i = 0; i < len; ++i)
    {
        array1[i] = new double[len];
    }

    transponse(len, second);
    for (int i = 0; i < len; ++i)
    {
        for (int j = 0; j < len; ++j)
        {
            array1[i][j] = mult_lines(len, first[i], second[j]);
        }
    }

    for (int i = 0; i < len; ++i)
    {
        for (int j = 0; j < len; ++j)
        {
            first[i][j] = array1[i][j];
        }
    }

    transponse(len, second);
    delete[] array1;
}

void mult_matrix_col (int len, double **first, double *second) {
    double *help = new double[len];
    for (int i = 0; i < len; ++i)
    {
        help[i] = mult_lines(len, first[i], second);
    }
    for (int i = 0; i < len; ++i)
    {
        second[i] = help[i];
    }
    delete[] help;
}

state mult_matrix_state (double **first, state U) {
    state ans;
    ans.rho = first[0][0]*U.rho + first[0][1]*U.rho_u + first[0][2]*U.rho_e;
    ans.rho_u = first[1][0]*U.rho + first[1][1]*U.rho_u + first[1][2]*U.rho_e;
    ans.rho_e = first[2][0]*U.rho + first[2][1]*U.rho_u + first[2][2]*U.rho_e;
    return ans;
}

state operator+ (state u1, state u2) {
    state ans;
    ans.rho = u1.rho + u2.rho;
    ans.rho_u = u1.rho_u + u2.rho_u;
    ans.rho_e = u1.rho_e + u2.rho_e;
    return ans;
}

state operator- (state u1, state u2) {
    state ans;
    ans.rho = u1.rho - u2.rho;
    ans.rho_u = u1.rho_u - u2.rho_u;
    ans.rho_e = u1.rho_e - u2.rho_e;
    return ans;
}

state operator* (double var, state u) {
    state ans;
    ans.rho = var*u.rho;
    ans.rho_u = var*u.rho_u;
    ans.rho_e = var*u.rho_e;
    return ans;
}

void print_matrix (int str, int col, double **array1) {
    for (int i = 0; i < str; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            std::cout << array1[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    state UL, UR;
    UL.rho = rho_l;
    UL.rho_u = rho_l * v_l;
    UL.rho_e = p_l / (gamma - 1.0);

    UR.rho = rho_r;
    UR.rho_u = rho_r * v_r;
    UR.rho_e = p_r / (gamma - 1.0);

    //tau = h * CFL / std::max(fabs(find_c(UR)), fabs(find_c(UL)));
    int NT = 0;

    int n = 3;  //matrices for present state
    double** omega = new double*[n];
    for (int i = 0; i < n; ++i)
    {
        omega[i] = new double[n];
    }
    double** reverse_omega = new double*[n];
    for (int i = 0; i < n; ++i)
    {
        reverse_omega[i] = new double[n];
    }
    double** lambda = new double*[n];
    for (int i = 0; i < n; ++i)
    {
        lambda[i] = new double[n];
    }
    double** mod_lambda = new double*[n];
    for (int i = 0; i < n; ++i)
    {
        mod_lambda[i] = new double[n];
    }
    double** A_matrix = new double*[n];
    for (int i = 0; i < n; ++i)
    {
        A_matrix[i] = new double[n];
    }
    double** mod_A_matrix = new double*[n];
    for (int i = 0; i < n; ++i)
    {
        mod_A_matrix[i] = new double[n];
    }

    state *layer1 = new state[NX];
    state *layer2 = new state[NX];
    for (int i = 0; i < NX / 2; ++i)
    {
        layer1[i] = UL;
    }
    for (int i = NX / 2; i < NX; ++i)
    {
        layer1[i] = UR;
    }

    double t = 0.0;
    for (NT = 0; t < T; ++NT)
    {
        t += tau;
        double max_lambda = 0.0;
        for(int j = 0; j < NX; ++j)
        {
            max_lambda = std::max(max_lambda, find_c(layer1[j]) + layer1[j].rho_u/layer1[j].rho);
        }
        tau = h * CFL / max_lambda * 0.1;
        tau = std::min(tau, 0.00001);

        if (NT % 2 == 0)
        {
            for (int j = 1; j < NX-1; ++j)
            {
                omega_U(layer1[j], omega);
                //print_matrix(3, 3, omega);+
                reverse_omega_U(layer1[j], reverse_omega);
                lambda_U(layer1[j], lambda);
                mod_lambda_U(layer1[j], mod_lambda);
                //construct A: A = reverse_omega*lambda*omega
                reverse_omega_U(layer1[j], A_matrix);
                mult_matrix(n, A_matrix, lambda);
                mult_matrix(n, A_matrix, omega);
                //print_matrix(3, 3, reverse_omega);+
                reverse_omega_U(layer1[j], mod_A_matrix);
                mult_matrix(n, mod_A_matrix, mod_lambda);
                mult_matrix(n, mod_A_matrix, omega);
                //print_matrix(3, 3, lambda);
                layer2[j] = layer1[j] - (0.5*tau/h) * mult_matrix_state(A_matrix, (layer1[j+1] - layer1[j-1])) + (0.5*tau/h) * mult_matrix_state(mod_A_matrix, (layer1[j+1] - 2.0*layer1[j] + layer1[j-1]));
            }
            layer2[0] = layer2[1];
            layer2[NX-1] = layer2[NX-2];
            //std::cout << layer2->rho;
        }
        else
        {
            for (int j = 1; j < NX-1; ++j)
            {
                omega_U(layer2[j], omega);
                reverse_omega_U(layer2[j], reverse_omega);
                lambda_U(layer2[j], lambda);
                mod_lambda_U(layer2[j], mod_lambda);
                //construct A: A = reverse_omega*lambda*omega
                reverse_omega_U(layer2[j], A_matrix);
                mult_matrix(n, A_matrix, lambda);
                mult_matrix(n, A_matrix, omega);

                reverse_omega_U(layer2[j], mod_A_matrix);
                mult_matrix(n, mod_A_matrix, mod_lambda);
                mult_matrix(n, mod_A_matrix, omega);

                layer1[j] = layer2[j] - (0.5*tau/h) * mult_matrix_state(A_matrix, (layer2[j+1] - layer2[j-1])) + (0.5*tau/h) * mult_matrix_state(mod_A_matrix, (layer2[j+1] - 2.0*layer2[j] + layer2[j-1]));

                if (t >= T / 2 and t <= T / 2 + tau and j == 50)
                {
                    print_matrix(3, 3, lambda);
                    print_matrix(3, 3, mod_lambda);
                    print_matrix(3, 3, omega);
                    print_matrix(3, 3, reverse_omega);
                    std::cout << layer1[50].rho << std::endl;
                    std::cout << layer1[50].rho_e << std::endl;
                    std::cout << layer1[50].rho_u << std::endl;
                }
            }
            layer1[0] = layer1[1];
            layer1[NX-1] = layer1[NX-2];
        }
        
    }
    std::cout << NT;
    //std::ofstream outfile;
    //outfile.open("CIR_scheme_rho.txt", std::ios_base::app);
    //if (NT % 2 == 0)
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << layer1[i].rho << " ";
    //    }
    //    outfile << layer1[NX-1].rho << std::endl;
    //}
    //else
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << layer2[i].rho << " ";
    //    }
    //    outfile << layer2[NX-1].rho << std::endl;
    //}
    //outfile.close();
    //std::cout << "rho was written" << std::endl;

    //outfile.open("CIR_scheme_u.txt", std::ios_base::app);
    //if (NT % 2 == 0)
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << layer1[i].rho_u / layer1[i].rho << " ";
    //    }
    //    outfile << layer1[NX-1].rho_u / layer1[NX-1].rho << std::endl;
    //}
    //else
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << layer2[i].rho_u / layer2[i].rho << " ";
    //    }
    //    outfile << layer2[NX-1].rho_u / layer2[NX-1].rho << std::endl;
    //}
    //outfile.close();
    //std::cout << "u was written" << std::endl;

    //outfile.open("CIR_scheme_e.txt", std::ios_base::app);
    //if (NT % 2 == 0)
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << layer1[i].rho_e / layer1[i].rho << " ";
    //    }
    //    outfile << layer1[NX-1].rho_e / layer1[NX-1].rho << std::endl;
    //}
    //else
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << layer2[i].rho_e / layer2[i].rho << " ";
    //    }
    //    outfile << layer2[NX-1].rho_e / layer2[NX-1].rho << std::endl;
    //}
    //outfile.close();
    //std::cout << "e was written" << std::endl;

    //outfile.open("CIR_scheme_p.txt", std::ios_base::app);
    //if (NT % 2 == 0)
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << pressure(layer1[i]) << " ";
    //    }
    //    outfile << pressure(layer1[NX-1]) << std::endl;
    //}
    //else
    //{
    //    for (int i = 0; i < NX - 1; ++i)
    //    {
    //        outfile << pressure(layer2[i]) << " ";
    //    }
    //    outfile << pressure(layer2[NX-1]) << std::endl;
    //}
    //outfile.close();
    //std::cout << "P was written" << std::endl;

    delete[] omega;
    delete[] reverse_omega;
    delete[] lambda;
    delete[] mod_lambda;
    delete[] A_matrix;
    delete[] mod_A_matrix;
return 0;
}
