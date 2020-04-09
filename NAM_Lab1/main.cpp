#include <iostream>
#include <cmath>
#include <iomanip>

double my_func(double x) {
    return sqrt(x + 2) - 2 * cos(x);
}

void dichotomy(double left, double right) {
    int i = 1;
    double center;
    double left_func;
    double right_func;
    double center_func;

    while(std::abs(right - left) > 0.01) {
        center = (left + right) / 2;
        left_func = my_func(left);
        right_func = my_func(right);
        center_func = my_func(center);

        std::cout << i << '\t' << std::fixed << std::setprecision(3) << left << '\t' << right << '\t' << center << '\t' << left_func
                  << '\t' << right_func << '\t' << center_func << '\t' << right - left << std::endl;

        i++;

        if (left_func * center_func < 0) {
            right = center;
        }
        else {
            left = center;
        }
    }

    center = (left + right) / 2;
    left_func = my_func(left);
    right_func = my_func(right);
    center_func = my_func(center);

    std::cout << i << '\t' << std::fixed << std::setprecision(3) << left << '\t' << right << '\t' << center << '\t' << left_func
              << '\t' << right_func << '\t' << center_func << '\t' << right - left << std::endl;
}

double sim_func(double x) {
    return x - (sqrt(x + 2) - 2 * cos(x)) / 2;
};

void simple_iteration_method(double x) {
    int i = 1;
    double new_x;
    double accuracy = 1;

    while (accuracy > 0.0000001) {
        new_x = sim_func(x);
        accuracy = std::abs(new_x - x);

        std::cout << i << '\t' << std::fixed << std::setprecision(8) << new_x << '\t' << accuracy << std::endl;

        x = new_x;
        i++;
    }
}

double derivative_my_func(double x) {
    return 1. / (2 * sqrt(x + 2)) + 2 * sin(x);
}

void newton_method(double x) {
    int i = 1;
    double new_x;
    double accuracy = 1;

    while (accuracy > 0.0000001) {
        new_x = x - my_func(x) / derivative_my_func(x);
        accuracy = std::abs(new_x - x);

        std::cout << i << '\t' << std::fixed << std::setprecision(12) << new_x << '\t' << accuracy << std::endl;

        x = new_x;
        i++;
    }
}

void secant_method(double x, double old_x) {
    int i = 1;
    double new_x;
    double accuracy = 1;

    while (accuracy > 0.0000001) {
        new_x = x - my_func(x) * (x - old_x) / (my_func(x) - my_func(old_x));
        accuracy = std::abs(new_x - x);

        std::cout << i << '\t' << std::fixed << std::setprecision(12) << new_x << '\t' << accuracy << std::endl;

        old_x = x;
        x = new_x;
        i++;
    }
}

int main() {
    dichotomy(0, 1);

    std::cout << std::endl;

    simple_iteration_method(0.63);

    std::cout << std::endl;

    newton_method(0.63);

    std:: cout << std::endl;

    secant_method(0.63, 0.62716377);

    return 0;
}
