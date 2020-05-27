#include <iostream>
#include <cmath>
#include <iomanip>

struct Polynomial {
    int dim;
    double* odds;

    explicit Polynomial(int dim, double odd1 = 0, double odd2 = 0) {
        this->dim = dim;
        this->odds = new double[this->dim + 1];
        this->odds[0] = odd1;

        if (this->dim == 1) {
            this->odds[1] = odd2;
        }
    }

    Polynomial operator+(const Polynomial& other) {
        if (this->dim >= other.dim) {
            Polynomial result(this->dim);

            for (int i = 0; i <= other.dim; i++) {
                result.odds[i] = this->odds[i] + other.odds[i];
            }

            for (int i = other.dim + 1; i <= this->dim; i++) {
                result.odds[i] = this->odds[i];
            }

            return result;
        } else {
            Polynomial result(other.dim);

            for (int i = 0; i <= this->dim; i++) {
                result.odds[i] = this->odds[i] + other.odds[i];
            }

            for (int i = this->dim + 1; i <= other.dim; i++) {
                result.odds[i] = other.odds[i];
            }

            return result;
        }
    }

    Polynomial operator*(const Polynomial& other) {
        Polynomial result(this->dim + other.dim);

        for (int i = 0; i <= this->dim; i++) {
            for (int j = 0; j <= other.dim; j++) {
                result.odds[i + j] += this->odds[i] * other.odds[j];
            }
        }

        return result;
    }

    friend std::ostream& operator<<(std::ostream& os, Polynomial& polynomial);
};

std::ostream& operator<<(std::ostream& os, Polynomial& polynomial) {
    for (int i = polynomial.dim; i > 0; i--) {
        os << std::fixed << std::setprecision(8) << polynomial.odds[i] << "x^{" << i << "} + ";
    }

    os << std::fixed << std::setprecision(8) << polynomial.odds[0];

    return os;
}

const int START = -2;
const int END = 2;

void set_uniform_nodes(double* nodes, int number) {
    for (int i = 0; i <= number; i++) {
        nodes[i] = START + 1. * i * (END - START) / number;
    }
}

void set_optimal_nodes(double* nodes, int number) {
    for (int i = 0; i <= number; i++) {
        nodes[i] = 1. * (START + END) / 2 + 1. * (END - START) / 2 * cos((2 * i + 1) * M_PI / (2 * (number + 1)));
    }
}

double func1(double x) {
    return log(x * x + 1) * sin(2 * x);
}

double func2(double x) {
    return sin(std::abs(x)) + 1;
}

Polynomial get_lagrange(int number, const double* nodes, int index) {
    Polynomial result(0, 1);

    for (int i = 0; i <= number; i++) {
        if (i != index) {
            Polynomial cur(1, (nodes[i] / (nodes[i] - nodes[index])), 1. / (nodes[index] - nodes[i]));

            result = result * cur;
        }
    }

    return result;
}

Polynomial get(int number, int num_func, const double* nodes) {
    Polynomial result(0);

    if (num_func == 1) {
        for (int i = 0; i <= number; i++) {
            result = result + Polynomial(0, func1(nodes[i])) * get_lagrange(number, nodes, i);
        }
    } else {
        for (int i = 0; i <= number; i++) {
            result = result + Polynomial(0, func2(nodes[i])) * get_lagrange(number, nodes, i);
        }
    }

    return result;
}

int main() {
    auto* nodes = new double[21];
    int dim;
    Polynomial polynomial_func1(0);
    Polynomial polynomial_func2(0);

    std::cout << "Optimal nodes:" << std::endl;

    dim = 3;
    set_optimal_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 5;
    set_optimal_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 7;
    set_optimal_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 10;
    set_optimal_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 15;
    set_optimal_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 20;
    set_optimal_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    std::cout << "Evenly spaced nodes:" << std::endl;

    dim = 3;
    set_uniform_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 5;
    set_uniform_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 7;
    set_uniform_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 10;
    set_uniform_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 15;
    set_uniform_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    dim = 20;
    set_uniform_nodes(nodes, dim);
    polynomial_func1 = get(dim, 1, nodes);
    polynomial_func2 = get(dim, 2, nodes);
    std::cout << dim << std::endl << polynomial_func1 << std::endl << polynomial_func2 << std::endl << std::endl;

    return 0;
}
