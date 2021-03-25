#include <iostream>
#include "utils.h"
#include <gmpxx.h>
#include <math.h>
#include "quad_sieve.h"
#include <vector>

using namespace std;

void quad_sieve(const mpz_class &n, mpz_class &fact1, mpz_class &fact2) {
    mpz_class smooth_bound = 19;
    long long primes_lt_bound = primes_lt_b(smooth_bound);

    mpz_class a[primes_lt_bound];
    mpz_class b[primes_lt_bound];

    mpz_class sqrtn = sqrt(n);

    for (long long i = 0; i <= primes_lt_bound; ++i) {
//        a[i] = i + sqrtn;
//        b[i] = a[i]*a[i] - n;
//
//        if ()
    }
}

/// REF https://www.geeksforgeeks.org/p-smooth-numbers-p-friable-number/
bool is_b_smooth(mpz_class n, mpz_class &b) {
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
    return  maximum <= b;
}


/// REF https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes
long long primes_lt_b(mpz_class &b) {
    long long count = 0;
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

    for(auto e : primes) {
        if (e) ++count;
    }
    return count;
}