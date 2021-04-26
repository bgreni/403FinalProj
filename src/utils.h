#pragma once

#include <iostream>
#include <gmpxx.h>
#include <gmp.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <iterator>
#include <map>
#include <utility>

#define Vec vector<mpz_class>
#define NOW chrono::high_resolution_clock::now()
#define DUR(x,y) chrono::duration_cast<chrono::microseconds>(y - x).count();

/// FOR COLORING OUTPUT TEXT
/// Just so I have this https://stackoverflow.com/questions/9158150/colored-output-in-c/9158263
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"

using namespace std;
using ull = unsigned long long;
using Matrix = vector<Vec>;
using SolRows = vector<pair<Vec, size_t>>;

inline mpz_class prod(Vec &vals) {
    mpz_class ret = vals[0];
    for (auto it = vals.begin()+1; it != vals.end(); ++it)
        ret *= *it;
    return ret;
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
        if (j == (long)exp_counts.size()) break;

        if (exp_counts[j] % 2 == 0)
            indices.push_back(j);
    }
    return indices;
}

void tonelli_shanks(const mpz_class n, mpz_class &p, mpz_class &x, mpz_class &other) {

    mpz_class q = p - 1;
    mpz_class s = 0;

    while (q % 2 == 0) {
        
    }

}

map<mpz_class, int> get_p_factors(mpz_class n, Vec base) {
    map<mpz_class, int> factors;
    int count = 0;

    if (n < 0)
        factors[-1] = 1;
    
    for (auto prime : base) {
        if (prime == -1) continue;

        while (n % prime == 0 && n != 0) {
            ++count;
            n /= prime;
        }
        factors[prime] = count;
        count = 0;
    }

    return factors;
}

Matrix transpose(Matrix &m) {
    Matrix new_matrix;

    for (size_t i = 0; i < m[0].size(); ++i) {
        Vec new_row;
        for (auto row : m) {
            new_row.push_back(row[i]);
        }
        new_matrix.push_back(new_row);
    }
    return new_matrix;
}

template<typename T>
void show(const T& t) {
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, ", "));
    cout << endl;
}