#include <iostream>
#include <cmath>
#include <iomanip>

double sim_func1(double x1, double x2) {
    return x1 - 1. * 500 / 891 * (x1 + x2 * cos(x1 - 2) - 0.5) + 1. * 10 / 891 * (x1 * x1 + x2 - 1);
}

double sim_func2(double x1, double x2) {
    return x2 + 1. * 50 / 99 * (x1 + x2 * cos(x1 - 2) - 0.5) - 1. * 100 / 99 * (x1 * x1 + x2 - 1);
}

void simple_iteration_method(double x1, double x2) {
    int i = 1;
    double new_x1;
    double new_x2;
    double accuracy = 1;

    while (accuracy >= 0.00001) {
        new_x1 = sim_func1(x1, x2);
        new_x2 = sim_func2(x1, x2);
        accuracy = std::max(std::abs(new_x1 - x1), std::abs(new_x2 - x2));

        std::cout << i << '\t' << std::fixed << std::setprecision(8) << new_x1 << '\t' << new_x2 << '\t' << accuracy << std::endl;

        x1 = new_x1;
        x2 = new_x2;
        i++;
    }
}

void seidels_method(double x1, double x2) {
    int i = 1;
    double new_x1;
    double new_x2;
    double accuracy = 1;

    while (accuracy >= 0.00001) {
        new_x1 = sim_func1(x1, x2);
        new_x2 = sim_func2(new_x1, x2);
        accuracy = std::max(std::abs(new_x1 - x1), std::abs(new_x2 - x2));

        std::cout << i << '\t' << std::fixed << std::setprecision(8) << new_x1 << '\t' << new_x2 << '\t' << accuracy << std::endl;

        x1 = new_x1;
        x2 = new_x2;
        i++;
    }
}

double nm_dfunc1(double x1, double x2) {
    return 1 - x2 * sin(x1 - 2);
}

double nm_dfunc2(double x1, double x2) {
    return cos(x1 - 2);
}

double nm_dfunc3(double x1, double x2) {
    return 2 * x1;
}

double nm_dfunc4(double x1, double x2) {
    return 1;
}

double nm_func1(double x1, double x2) {
    return x1 + x2 * cos(x1 - 2) - 0.5;
}

double nm_func2(double x1, double x2) {
    return x1 * x1 + x2 - 1;
}

void gauss_method(double** matrix, double* f, double* result, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            matrix[i][j] /= matrix[i][i];
        }

        f[i] /= matrix[i][i];

        for (int j = i + 1; j < n; j++) {
            double c = -1 * matrix[j][i];
            matrix[j][i] = 0;

            for (int k = i + 1; k < n; k++) {
                matrix[j][k] += matrix[i][k] * c;
            }
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        result[i] = f[i];

        for (int j = n - 1; j > i; j--) {
            result[i] -= matrix[i][j] * result[j];
        }
    }
}

void newtons_method_for_matrix(double x1, double x2) {
    auto** jacobi_matrix = new double*[2];
    auto* f = new double[2];
    int i = 1;
    double new_x1;
    double new_x2;
    double accuracy = 1;

    for (int j = 0; j < 2; j++) {
        jacobi_matrix[j] = new double[2];
    }

    while (accuracy >= 0.00001) {
        jacobi_matrix[0][0] = nm_dfunc1(x1, x2);
        jacobi_matrix[0][1] = nm_dfunc2(x1, x2);
        jacobi_matrix[1][0] = nm_dfunc3(x1, x2);
        jacobi_matrix[1][1] = nm_dfunc4(x1, x2);
        f[0] = -1 * nm_func1(x1, x2);
        f[1] = -1 * nm_func2(x1, x2);

        auto* result = new double[2];

        gauss_method(jacobi_matrix, f, result, 2);

        new_x1 = result[0] + x1;
        new_x2 = result[1] + x2;

        accuracy = std::max(std::abs(new_x1 - x1), std::abs(new_x2 - x2));

        std::cout << i << '\t' << std::fixed << std::setprecision(8) << new_x1 << '\t' << new_x2 << '\t' << accuracy << std::endl;

        x1 = new_x1;
        x2 = new_x2;
        i++;
    }
}

double my_func1(double x, double c) {
    return x + c * cos(x - 2) - 0.5;
}

double derivative_my_func1(double x, double c) {
    return 1 - c * sin(x - 2);
}

double my_func2(double x, double c) {
    return c * c + x - 1;
}

double derivative_my_func2(double x, double c) {
    return 1;
}

double newtons_method(double x, double c, int num) {
    double new_x = 0;
    double accuracy = 1;

    while (accuracy > 0.00001) {
        switch (num) {
            case 1:
                new_x = x - my_func1(x, c) / derivative_my_func1(x, c);
                break;
            case 2:
                new_x = x - my_func2(x, c) / derivative_my_func2(x, c);
                break;
        }

        accuracy = std::abs(new_x - x);
        x = new_x;
    }

    return x;
}

void gauss_seidels_method(double x1, double x2) {
    int i = 1;
    double new_x1;
    double new_x2;
    double accuracy = 1;

    while (accuracy > 0.00001) {
        new_x1 = newtons_method(x1, x2, 1);
        new_x2 = newtons_method(x2, new_x1, 2);

        accuracy = std::max(std::abs(new_x1 - x1), std::abs(new_x2 - x2));

        std::cout << i << '\t' << std::fixed << std::setprecision(8) << new_x1 << '\t' << new_x2 << '\t' << accuracy << std::endl;

        x1 = new_x1;
        x2 = new_x2;
        i++;
    }
}

int main() {
    simple_iteration_method(0.45, 0.8);

    std::cout << std::endl;

    seidels_method(0.45, 0.8);

    std::cout << std::endl;

    newtons_method_for_matrix(0.45, 0.8);

    std::cout << std::endl;

    gauss_seidels_method(0.45, 0.8);

    std::cout << std::endl;

    return 0;
}
