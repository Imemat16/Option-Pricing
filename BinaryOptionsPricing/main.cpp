#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

int index(int t, int i) {
    return t * (t + 1) / 2 + i;
}

double normalCDF(double x) {
    return 0.5 * erfc(-x * M_SQRT1_2); // Standard normal CDF using error function
}

double blackScholesPrice(double S, double K, double T, double r, double sigma, bool isCall) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    if (isCall)
        return S * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
    else
        return K * exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1);
}

void exportTreeCSV(const vector<double>& tree, int N, const string& filename) {
    ofstream file(filename);
    for (int t = 0; t <= N; ++t) {
        for (int i = 0; i <= t; ++i) {
            file << tree[index(t, i)];
            if (i != t) file << ",";
        }
        file << "\n";
    }
    file.close();
}

void printTree(const vector<double>& tree, int N, const string& label) {
    cout << "\n" << label << " Tree:\n";
    for (int t = 0; t <= N; ++t) {
        cout << "t=" << t << ": ";
        for (int i = 0; i <= t; ++i) {
            cout << fixed << setprecision(2) << tree[index(t, i)] << " ";
        }
        cout << "\n";
    }
}

double priceOptionBinomial(
    double S0, double K, double T, double r, double sigma,
    int N, bool isCall, bool isAmerican
    ) {
    double dt = T / N;
    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p = (exp(r * dt) - d) / (u - d);

    int totalNodes = (N + 1) * (N + 2) / 2;
    vector<double> stockTree(totalNodes);
    vector<double> optionTree(totalNodes);

    // Fill stock price tree
    for (int t = 0; t <= N; ++t) {
        for (int i = 0; i <= t; ++i) {
            stockTree[index(t, i)] = S0 * pow(u, i) * pow(d, t - i);
        }
    }

    // Terminal payoffs
    for (int i = 0; i <= N; ++i) {
        double ST = stockTree[index(N, i)];
        optionTree[index(N, i)] = isCall ? max(0.0, ST - K) : max(0.0, K - ST);
    }

    // Backward induction
    for (int t = N - 1; t >= 0; --t) {
        for (int i = 0; i <= t; ++i) {
            double hold = exp(-r * dt) * (p * optionTree[index(t + 1, i + 1)] + (1 - p) * optionTree[index(t + 1, i)]);
            if (isAmerican) {
                double exercise = isCall ? max(0.0, stockTree[index(t, i)] - K) : max(0.0, K - stockTree[index(t, i)]);
                optionTree[index(t, i)] = max(hold, exercise);
            } else {
                optionTree[index(t, i)] = hold;
            }
        }
    }

    // Optional: print and export trees
    printTree(stockTree, N, "Stock Price");
    printTree(optionTree, N, "Option Price");
    exportTreeCSV(stockTree, N, "stock_tree.csv");
    exportTreeCSV(optionTree, N, "option_tree.csv");

    return optionTree[0];
}

int main() {
    // Input parameters
    double S0 = 100.0;
    double K = 100.0;
    double T = 1.0;
    double r = 0.05;
    double sigma = 0.2;
    int N = 50;
    bool isCall = true;
    bool isAmerican = false;

    // Price using binomial tree
    double price = priceOptionBinomial(S0, K, T, r, sigma, N, isCall, isAmerican);

    // Compare with Black-Scholes
    double bsPrice = blackScholesPrice(S0, K, T, r, sigma, isCall);

    // Output results
    cout << "\nOption Type: " << (isCall ? "Call" : "Put") << " | "
         << (isAmerican ? "American" : "European") << "\n";
    cout << "Binomial Tree Price: " << fixed << setprecision(4) << price << "\n";
    cout << "Black-Scholes Price: " << fixed << setprecision(4) << bsPrice << "\n";

    return 0;
}
