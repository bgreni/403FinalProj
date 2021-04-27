#include <gmpxx.h>
#include <iostream>
#include "quad_sieve.h"
#include "utils.h"
#include <assert.h>
#include "test_numbers.h"

using namespace std;

int main() {
    mpz_class n = RSA_90bit;
    mpz_class f1, f2;
//    mpz_pow_ui(n.get_mpz_t(), n.get_mpz_t(), 10);
//    n = sqrt(n);
//    double m = 1000;
//    m = sqrt(m);
//    get_factor_base(n,n);
    QSFact qs = QSFact();
    auto s = NOW;
    qs.quad_sieve(n, f1, f2);
    auto e = NOW;
    auto t = DUR(s, e);
    mpz_class res = f1*f2;
    ASSERT((res == n), n, res);
    cout << "time taken (seconds): " << t << endl;
    cout << "factor 1: " << f1 << endl;
    cout << "factor 2: " << f2 << endl;
}