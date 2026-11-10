#include "pricing.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace options {
namespace {

void validate(const Contract& c) {
    if (!std::isfinite(c.spot) || !std::isfinite(c.strike) ||
        !std::isfinite(c.years) || !std::isfinite(c.rate) ||
        !std::isfinite(c.volatility) || c.spot <= 0 || c.strike <= 0 ||
        c.years < 0 || c.volatility < 0) {
        throw std::invalid_argument("Use positive spot and strike, nonnegative expiry and volatility, and finite inputs.");
    }
}

double intrinsic(const Contract& c, double spot) {
    const double spread = spot - c.strike;
    return std::max(0.0, c.payoff == Payoff::call ? spread : -spread);
}

double normal_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

}  // namespace

double black_scholes(const Contract& c) {
    validate(c);
    if (c.years == 0) {
        return intrinsic(c, c.spot);
    }
    const double discounted_strike = c.strike * std::exp(-c.rate * c.years);
    if (c.volatility == 0) {
        const double spread = c.spot - discounted_strike;
        return std::max(0.0, c.payoff == Payoff::call ? spread : -spread);
    }
    const double scale = c.volatility * std::sqrt(c.years);
    const double d1 = (std::log(c.spot) - std::log(c.strike) + c.rate * c.years) / scale + scale / 2;
    const double d2 = d1 - scale;
    if (c.payoff == Payoff::call) {
        return c.spot * normal_cdf(d1) - discounted_strike * normal_cdf(d2);
    }
    return discounted_strike * normal_cdf(-d2) - c.spot * normal_cdf(-d1);
}

double binomial(const Contract& c, int steps, Exercise exercise, Lattice* lattice) {
    validate(c);
    if (steps <= 0) {
        throw std::invalid_argument("Steps must be positive.");
    }
    if (lattice && steps > 2000) {
        throw std::invalid_argument("Tree export supports at most 2000 steps.");
    }
    const double dt = c.years / steps;
    const double jump = c.volatility * std::sqrt(dt);
    const double drift = c.rate * dt;
    const double discount = std::exp(-drift);
    const bool deterministic = jump == 0.0;
    // CRR uses reciprocal up and down factors. Shift both paths by the
    // risk-free drift only in the zero-volatility limit.
    const double center = deterministic ? drift : 0.0;
    const double p = deterministic ? 1.0 : std::expm1(drift + jump) / std::expm1(2 * jump);
    if (!std::isfinite(p) || p < 0 || p > 1 || !std::isfinite(discount)) {
        throw std::domain_error("Invalid risk-neutral tree. Try more steps or less extreme inputs.");
    }

    const double log_spot = std::log(c.spot);
    auto spot_at = [&](int level, int up_moves) {
        const double spot = std::exp(log_spot + level * center + (2.0 * up_moves - level) * jump);
        if (!std::isfinite(spot)) {
            throw std::overflow_error("A stock price in the tree exceeds the numeric range.");
        }
        return spot;
    };

    const auto level_count = static_cast<std::size_t>(steps) + 1;
    std::vector<double> values(level_count);
    if (lattice) {
        lattice->stock.resize(level_count);
        lattice->option.resize(level_count);
    }
    auto record = [&](int level) {
        if (lattice) {
            lattice->option[level].assign(values.begin(), values.begin() + level + 1);
            auto& stock = lattice->stock[level];
            stock.resize(static_cast<std::size_t>(level) + 1);
            for (int j = 0; j <= level; ++j) {
                stock[j] = spot_at(level, j);
            }
        }
    };

    for (int j = 0; j < steps; ++j) {
        values[j] = intrinsic(c, spot_at(steps, j));
    }
    values.back() = intrinsic(c, spot_at(steps, steps));
    record(steps);
    for (int level = steps - 1; level >= 0; --level) {
        // Moving left to right leaves the next value at j + 1 untouched.
        for (int j = 0; j <= level; ++j) {
            const double continuation = discount * ((1 - p) * values[j] + p * values[j + 1]);
            values[j] = exercise == Exercise::american
                ? std::max(continuation, intrinsic(c, spot_at(level, j)))
                : continuation;
        }
        record(level);
    }
    if (!std::isfinite(values.front())) {
        throw std::overflow_error("The option price exceeds the numeric range.");
    }
    return values.front();
}

}  // namespace options
