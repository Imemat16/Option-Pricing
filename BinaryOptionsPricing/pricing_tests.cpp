#include "pricing.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

}  // namespace

int main() {
    using namespace options;
    try {
        const Contract call{100, 100, 1, 0.05, 0.20, Payoff::call};
        Contract put = call;
        put.payoff = Payoff::put;
        const double reference = black_scholes(call);
        require(std::abs(reference - 10.4505835722) < 1e-9, "Black-Scholes benchmark");
        require(std::abs(reference - black_scholes(put) - (100 - 100 * std::exp(-0.05))) < 1e-9,
                "European put-call parity");

        double last_error = 1.0;
        for (int steps : {50, 200, 800}) {
            const double tree = binomial(call, steps, Exercise::european);
            const double error = std::abs(tree - reference);
            require(error < last_error, "European call convergence on an even step grid");
            last_error = error;
        }
        require(last_error < 0.003, "800-step pricing tolerance");

        const double european_put = binomial(put, 200, Exercise::european);
        const double american_put = binomial(put, 200, Exercise::american);
        require(american_put > european_put, "American put early exercise premium");
        require(std::abs(binomial(call, 200, Exercise::american) -
                         binomial(call, 200, Exercise::european)) < 1e-9,
                "Non-dividend-paying American call matches European call");

        Contract expired = put;
        expired.spot = 90;
        expired.years = 0;
        require(black_scholes(expired) == 10, "Expired Black-Scholes put");
        require(std::abs(binomial(expired, 10, Exercise::american) - 10) < 1e-9,
                "Expired binomial put");

        Contract deterministic = call;
        deterministic.volatility = 0;
        require(std::abs(binomial(deterministic, 50, Exercise::european) -
                         black_scholes(deterministic)) < 1e-9, "Zero-volatility limit");

        Lattice lattice;
        const double exported = binomial(put, 4, Exercise::american, &lattice);
        require(lattice.stock.size() == 5 && lattice.stock.back().size() == 5 &&
                lattice.option.front().size() == 1 && lattice.option.front().front() == exported &&
                exported == binomial(put, 4, Exercise::american), "Export matches rolling price");

        bool rejected_steps = false;
        try {
            binomial(call, 0, Exercise::european);
        } catch (const std::invalid_argument&) {
            rejected_steps = true;
        }
        require(rejected_steps, "Reject zero steps");

        Contract invalid_tree = call;
        invalid_tree.rate = 1.0;
        bool rejected_probability = false;
        try {
            binomial(invalid_tree, 1, Exercise::european);
        } catch (const std::domain_error&) {
            rejected_probability = true;
        }
        require(rejected_probability, "Reject out-of-range risk-neutral probability");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
