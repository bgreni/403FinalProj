#pragma once

#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <iterator>

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

/// REF https://www.geeksforgeeks.org/p-smooth-numbers-p-friable-number/
inline bool is_b_smooth(mpz_class n, mpz_class &b) {
    mpz_class maximum = -1;

    while (mpz_even_p(n.get_mpz_t())) {
        if (2 > maximum)
            maximum = 2;
        n /= 2;
    }

    for (int i = 3; i < sqrt(n) ; i += 2) {
        while (n % i == 0) {
            if (i > maximum)
                maximum = i;
            n /= i;
        }
    }

    if (n > 2 && n > maximum)
        maximum = n;
//    cout << "MAX " << maximum << endl;
    return  maximum <= b;
}

inline vector<long> get_even_indices(const vector<long> &exp_counts, long primes_lt_count) {
    vector<long> indices;

    for (long j = 0; j <= primes_lt_count; ++j) {
//        if (j == (long)exp_counts.size()) break;

        if (exp_counts[j] % 2 == 0)
            indices.push_back(j);
    }
    return indices;
}


template<typename T>
void show(const T& t) {
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, ", "));
    cout << endl;
}