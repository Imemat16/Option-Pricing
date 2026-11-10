#include "pricing.hpp"

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

void usage(std::ostream& stream) {
    stream << "Usage: option-pricer [spot strike years rate volatility steps call|put european|american [--export-tree]]\n"
           << "Rates and volatility are decimals. Run without arguments for a sample.\n";
}

double number(const char* text) {
    std::size_t parsed = 0;
    const double result = std::stod(text, &parsed);
    if (text[parsed] != '\0') {
        throw std::invalid_argument("Invalid numeric argument.");
    }
    return result;
}

int step_count(const char* text) {
    std::size_t parsed = 0;
    const int result = std::stoi(text, &parsed);
    if (text[parsed] != '\0') {
        throw std::invalid_argument("Steps must be an integer.");
    }
    return result;
}

void write_csv(const char* filename, const std::vector<std::vector<double>>& levels) {
    std::ofstream out(filename);
    out << std::setprecision(17);
    for (const auto& level : levels) {
        out << level.front();
        for (std::size_t j = 1; j < level.size(); ++j) {
            out << ',' << level[j];
        }
        out << '\n';
    }
    out.close();
    if (!out) {
        throw std::runtime_error(std::string("Could not write ") + filename);
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    using namespace options;
    if (argc == 2 && std::string_view(argv[1]) == "--help") {
        usage(std::cout);
        return 0;
    }
    if (argc != 1 && argc != 9 && argc != 10) {
        usage(std::cerr);
        return 2;
    }
    try {
        Contract contract{100.0, 100.0, 1.0, 0.05, 0.20, Payoff::call};
        int steps = 50;
        Exercise exercise = Exercise::european;
        const bool export_tree = argc == 10;
        if (argc > 1) {
            const std::string_view payoff = argv[7];
            const std::string_view style = argv[8];
            if ((payoff != "call" && payoff != "put") ||
                (style != "european" && style != "american") ||
                (export_tree && std::string_view(argv[9]) != "--export-tree")) {
                throw std::invalid_argument("Use call or put, european or american, and optionally --export-tree.");
            }
            contract = {number(argv[1]), number(argv[2]), number(argv[3]), number(argv[4]),
                        number(argv[5]), payoff == "call" ? Payoff::call : Payoff::put};
            steps = step_count(argv[6]);
            exercise = style == "european" ? Exercise::european : Exercise::american;
        }

        Lattice lattice;
        const double price = binomial(contract, steps, exercise, export_tree ? &lattice : nullptr);
        std::cout << std::fixed << std::setprecision(6)
                  << "Binomial (" << steps << " steps): " << price << '\n';
        if (exercise == Exercise::european) {
            std::cout << "Black-Scholes: " << black_scholes(contract) << '\n';
        }
        if (export_tree) {
            write_csv("stock_tree.csv", lattice.stock);
            write_csv("option_tree.csv", lattice.option);
            std::cout << "Wrote stock_tree.csv and option_tree.csv\n";
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
