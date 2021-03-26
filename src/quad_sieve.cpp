#include <iostream>
#include "utils.h"
#include <gmpxx.h>
#include <math.h>
#include "quad_sieve.h"
#include <vector>

using namespace std;

void quad_sieve(const mpz_class &n, mpz_class &fact1, mpz_class &fact2) {
    mpz_class smooth_bound;
    get_factor_base(n, smooth_bound);
    if (19 > smooth_bound)
        smooth_bound = 19;
    cout << "SMOOTH BOUND " << smooth_bound << endl;
    vector<long> primes;
    long primes_lt_count = primes_lt_b(smooth_bound, primes);
    cout << "NUM PRIMES: " << primes_lt_count << endl;

    Vec a(primes_lt_count+1);
    Vec b(primes_lt_count+1);

    mpz_class interval_lower_bound = sqrt(n)+1;
//    mpz_class interval_upper = sqrt(2*n) - 1;

    vector<vector<long>> matrix;

    for (long i = 0; i <= primes_lt_count; ++i) {
        // pray nothing does wrong here
        a[i] = i + interval_lower_bound;
        cout << "Ai: " << a[i] << endl;
        b[i] = (a[i]*a[i]) - n;
        cout << "Bi: " << b[i] << endl;
        if (is_b_smooth(b[i], smooth_bound)) {
            vector<long> exps_counts = get_factor_vector(b[i]);
            cout << "EXP COUNTS: ";
            show(exps_counts);
            matrix.push_back(exps_counts);
            vector<long> indices = get_even_indices(exps_counts, primes_lt_count);

            if (indices.size() == 0) continue;

            cout << "I: " << i << endl;
            cout << "INDICES: ";
            show(indices);

            mpz_class x = 1;

            for (auto j : indices) {
                x  = (x * a[j]) % n;
            }

            mpz_class y = 1;
            mpz_class res;
            for (auto j : indices) {
                mpz_class ex = 0;
                for (auto k : indices) {
                    ex += matrix[k][j];
                }
                mpz_powm(res.get_mpz_t(), mpz_class(primes[j]).get_mpz_t(), ex.get_mpz_t(), n.get_mpz_t());
                y *= res;
            }
            cout << "X: " << x << endl;
            cout << "Y: " << y << endl;
            if (x != (abs(y) % n) && (x != y))
                continue;
            else if (x == (abs(y) % n)) {
                fact1 = gcd(x - y, n);
                fact2 = gcd(x + y, n);
                return;
            }
        }
    }
}


/// REF https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
long primes_lt_b(mpz_class &b, vector<long> &prime_set) {

    // estimate of number of primes less than B
    // https://primes.utm.edu/howmany.html
    // only want to reserve roughly as much memory as we need
    long num_p_estimate = b.get_d() / (log(b.get_d()) - 1);
    prime_set.reserve(num_p_estimate);

    long count = 0;
    auto B = b.get_ui();
    vector<bool> primes(B, true);
    primes[0] = primes[1] = false;

    for (long unsigned i = 2; i < sqrt(B); ++i) {
        if (primes[i]) {
            for (long unsigned j = i*i; j < B; j+=i) {
                primes[j] = false;
            }
        }
    }

//    int k = 0;
    for(unsigned long i = 0; i < primes.size(); ++i) {
        if (primes[i] && mpz_legendre(b.get_mpz_t(), mpz_class(primes[i]).get_mpz_t()) == 1) {
            ++count;
            prime_set.push_back(i);
        }
    }
    return count;
}

/// REF https://www.geeksforgeeks.org/print-all-prime-factors-of-a-given-number/
vector<long> get_factor_vector(mpz_class n) {
    /// TODO find some estimate of the number of elements we might need to store
    vector<long> factors_exps;
    long count = 0;
    cout << "FACTOR PRIMES: ";
    while (n % 2 == 0) {
        ++count;
        cout << 2 << " ";
        n /= 2;
    }
    if (count != 0)
        factors_exps.push_back(count);
    count = 0;

    // TODO not sure if using mpz as an iterator variable is a great idea
    for (mpz_class i = 3; i <= sqrt(n); i += 2) {
        while (n % i == 0) {
            ++count;
            n /= i;
            cout << i << " ";
        }
        if (n == 0) break;
        if (count != 0)
            factors_exps.push_back(count);
        count = 0;
    }

    if (n > 2) {
        factors_exps.push_back(n.get_si());
        cout << n << endl;
    }

    return factors_exps;
}