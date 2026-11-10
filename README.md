# Option Pricer

A C++17 command-line pricer for vanilla calls and puts on non-dividend-paying assets. It implements Black–Scholes for European exercise and a Cox–Ross–Rubinstein (CRR) binomial tree for European and American exercise. At each American node, the tree compares continuation value with immediate exercise.

## Build and test

Use CMake 3.16 or later and a C++17 compiler.

```sh
cmake -S BinaryOptionsPricing -B build
cmake --build build --config Release
cmake -E chdir build ctest -C Release --output-on-failure
```

## Run

With a single-configuration generator, run `./build/option-pricer`. With Visual Studio, run `.\build\Release\option-pricer.exe`. No arguments runs the 50-step European call example. To pass your own inputs, use this order:

```sh
./build/option-pricer 100 100 1 0.05 0.20 200 put american
./build/option-pricer 100 100 1 0.05 0.20 50 call european --export-tree
```

The inputs are spot, strike, years to expiration, annual continuously compounded risk-free rate, annual volatility, tree steps, payoff type, and exercise style. Rates and volatility are entered as decimals. The export flag writes `stock_tree.csv` and `option_tree.csv` to the working directory. Each row is a time level, with columns ordered by increasing up moves. Export is capped at 2,000 steps.

## Convergence example

European call with spot and strike of 100, one year to expiry, a 5% rate, and 20% volatility. The Black–Scholes price is 10.450584.

| CRR steps | Price | Absolute error |
| ---: | ---: | ---: |
| 50 | 10.410692 | 0.039892 |
| 200 | 10.440591 | 0.009992 |
| 800 | 10.448084 | 0.002499 |

This is an even-step comparison. CRR convergence need not be monotone across every step count. The tests check convergence on this grid, European put-call parity, American early exercise, expiry, zero volatility, and invalid tree inputs.

Pricing uses backward induction with O(N²) time and O(N) memory. CSV export retains both full lattices and uses O(N²) memory.
