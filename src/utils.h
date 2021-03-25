#pragma once

#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <chrono>
#include <cmath>

#define Vec vector<mpz_class>
#define NOW chrono::high_resolution_clock::now()
#define DUR(x,y) chrono::duration_cast<chrono::microseconds>(y - x).count();

/// FOR COLORING OUTPUT TEXT
/// Just so I have this https://stackoverflow.com/questions/9158150/colored-output-in-c/9158263
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"

using namespace std;

inline void prod(Vec &vals, mpz_class &ret) {
    ret = vals[0];
    for (auto it = vals.begin()+1; it != vals.end(); ++it)
        ret *= *it;
}

inline void sum(Vec &vals, mpz_class &ret) {
    ret = 0;
    for (auto it = vals.begin(); it != vals.end(); ++it)
        ret += *it;
}

/// Formula taken from page 4 of this paper
/// http://www.damianball.com/pdf/portfolio/quadratic-sieve.pdf
inline void get_factor_base(const mpz_class &n, mpz_class &base) {
    // Want to do the arithmetic using floats, but wish to store
    // the result as an integer
    // likely to lose a reasonable amount of precision here at larger values of n
    // but it likely won't be important
    double N = n.get_d();
    base = pow(exp(sqrt(log(N) * log(log(N)))), (sqrt(2.0)/4.0));
}