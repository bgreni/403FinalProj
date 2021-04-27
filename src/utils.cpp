#include "utils.h"
#include <assert.h> 
#include "wrappers.h"


// https://gmplib.org/list-archives/gmp-devel/2006-May/000633.html
void tonelli_shanks(const mpz_class n, const mpz_class &p, mpz_class &sol1, mpz_class &sol2) {

    mpz_class q, r, z, c, t, t2, b;
    ulong s=0, i, m;

    q = p - 1;
    while (q % 2 == 0) {
        ++s;
        q /= 2;
    }
    
    if (s == 1) {
        powm(sol1, n, (p+1) / 4, p);
        sol2 = p - sol1;
        return;
    }
    z = 2;
    while (p - 1 != legendre(z, p) && z < p)
        ++z;

    powm(c, z, q, p);
    powm(r, n, (q+1) / 2, p);
    powm(t, n, q, p);
    m = s;


    while ((t - 1) % p != 0) {
        t2 = (t * t) % p;
        i = 1;
        for (; i < m; ++i) {
            if ((t2 - 1) % p == 0) break;
            t2 = (t2 * t2) % p;
        }

        powm_ui(b, c, 1 << (m - i - 1), p);
        r = (r * b) % p;
        c = (b * b) % p;
        t = (t * c) % p;
        m = i;
    }
    sol1 = r;
    sol2 = p - r;
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
        if (count > 0) {
            factors[prime] = count;
            count = 0;
        }
    }

    return factors;
}

Matrix transpose(Matrix &m) {
    Matrix new_matrix(m[0].size(), Vec(m.size()));

    for (size_t i = 0; i < m[0].size(); ++i) {
        for (size_t j = 0; j < m.size(); ++j) {
            new_matrix[i][j] = m[j][i];
        }
    }
    return new_matrix;
}

bool matrix_eq(const Matrix &m1, const Matrix &m2) {
    if (m1.size() != m2.size()) return false;

    for (size_t i = 0; i < m1.size(); ++i) {
        if (m1[i].size() != m2[i].size()) return false;
        for (size_t j = 0; j < m1[i].size(); ++j) {
            if (m1[i][j]  != m2[i][j]) return false;
        }
    }
    return true;
}