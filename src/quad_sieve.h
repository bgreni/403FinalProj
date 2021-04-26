#include <gmpxx.h>
#include <utils.h>


#pragma once

using namespace std;

class QSFact {
    mpz_class curr_x = 0;
    mpz_class root;
    long smooth_bound;
    mpz_class interval_size;
    Vec factor_base;

    public:
        void quad_sieve(const mpz_class &n, mpz_class &fact1, mpz_class &fact2);
        void initialize(const mpz_class &n);
        void create_factor_base(const mpz_class &n);
        void gen_smooth_numbers(const mpz_class &n, Vec &smooths, Vec xlist);
        Matrix gen_matrix(const Vec &smooths, const mpz_class &n);\
        void solve_linear(Matrix &matrix, vector<bool> &flagged, SolRows &solution_rows);
        void gauss(Matrix &A, vector<bool> &flagged);
        Vec find_dependencies(pair<Vec, size_t> &solution, Matrix &matrix, vector<bool> flagged);
};

