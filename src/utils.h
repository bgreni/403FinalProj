
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
#include <algorithm>
#include <sstream>

using namespace std;

#define Vec vector<mpz_class>
#define NOW chrono::high_resolution_clock::now()
#define DUR(x,y) chrono::duration<double, chrono::seconds::period>(y - x).count()

// test macros
#define ASSERT(x, y, z) { if (!x) {cout << RED << __FUNCTION__ << " failed on line " << __LINE__ << " "; \
    cout << "Expected " << y << " got " << z << RESET << endl; return false;} }
#define IS_TRUE(x) { if (!x) {cout << RED << __FUNCTION__ << " failed on line " << __LINE__ << RESET << endl; return false;} }
#define PASSED auto t = DUR(s,e); cout << GREEN << __FUNCTION__ << " passed in " << t << " seconds" << RESET << endl; return true

/// FOR COLORING OUTPUT TEXT
/// Just so I have this https://stackoverflow.com/questions/9158150/colored-output-in-c/9158263
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"

using ull = unsigned long long;
using Matrix = vector<Vec>;
using SolRows = vector<pair<Vec, size_t>>;


inline mpz_class prod(Vec &vals) {
    mpz_class ret = 1;
    for (auto it = vals.begin(); it != vals.end(); ++it)
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
// inline void get_factor_base(const mpz_class &n, mpz_class &base) {
//     // Want to do the arithmetic using floats, but wish to store
//     // the result as an integer
//     // likely to lose a reasonable amount of precision here at larger values of n
//     // but it likely won't be important
//     double N = n.get_d();
//     base = pow(exp(sqrt(log(N) * log(log(N)))), (sqrt(2.0)/4.0));
// }


template<typename T>
string mat_to_string(const vector<vector<T>>& t) {
    stringstream ss;
    for (auto row : t) {
        ss << "{ ";
        for (auto el : row) {
            ss << el << ", ";
        }
        ss << "}\n";
    }
    ss << endl;
    return ss.str();
}

template<typename T>
string vec_to_string(const vector<T>& t) {
    stringstream ss;
    ss << "{ ";
    for (auto el : t) {
        ss << el << ", ";
    }
    ss << "}\n";
    return ss.str();
}

template<typename T>
void show_mat(const vector<vector<T>>& t) {
    cout << mat_to_string(t) << endl;
}

template<typename T>
void showvec(const vector<T>& t) {
    cout << vec_to_string(t) << endl;
}

template<typename T>
bool vec_equality(const vector<T> &t1, const vector<T> &t2) {
    if (t1.size() != t2.size()) return false;
    for (size_t i = 0; i < t1.size(); ++i) {
        if (t1[i] != t2[i]) return false;
    }
    return true;
}

void tonelli_shanks(const mpz_class n, const mpz_class &p, mpz_class &x, mpz_class &other);
map<mpz_class, int> get_p_factors(mpz_class n, Vec base);
Matrix transpose(Matrix &m);
bool matrix_eq(const Matrix &m1, const Matrix &m2);
