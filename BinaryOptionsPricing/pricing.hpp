#pragma once

#include <vector>

namespace options {

enum class Payoff { call, put };
enum class Exercise { european, american };

struct Contract {
    double spot;
    double strike;
    double years;
    double rate;
    double volatility;
    Payoff payoff;
};

struct Lattice {
    std::vector<std::vector<double>> stock;
    std::vector<std::vector<double>> option;
};

double black_scholes(const Contract& contract);
double binomial(const Contract& contract, int steps, Exercise exercise, Lattice* lattice = nullptr);

}  // namespace options
