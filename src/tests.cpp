#include <iostream>
#include <gmpxx.h>
#include "utils.h"
#include <vector>
#include "quad_sieve.h"
#include <math.h>

using namespace std;
#define ASSERT(x, y, z) { if (!x) {cout << RED << __FUNCTION__ << " failed on line " << __LINE__ << " "; \
    cout << "Expected " << y << " got " << z << RESET << endl; return false;} }
#define IS_TRUE(x) { if (!x) {cout << RED << __FUNCTION__ << " failed on line " << __LINE__ << RESET << endl; return false;} }
#define PASSED cout << GREEN << __FUNCTION__ << " passed in " << t << " ms" << RESET << endl;

bool test_prod1() {
    auto s = NOW;

    mpz_class ret = 0;
    Vec vals = {1,2,3};
    prod(vals, ret);
    ASSERT((ret == 6), 6, ret);

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_prod2() {
    auto s = NOW;

    mpz_class ret = 0;
    Vec vals = {300,232,454};
    prod(vals, ret);
    ASSERT((ret == 31598400), 31598400, ret);

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_prod3() {
    auto s = NOW;
    mpz_class n = 10;
    mpz_pow_ui(n.get_mpz_t(), n.get_mpz_t(), 100);

    mpz_class ret = 0;
    Vec vals(100, 10);
    prod(vals, ret);
    ASSERT((ret == n), n, ret);

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_sum1() {
    auto s = NOW;

    mpz_class ret = 0;
    Vec vals = {1,2,3};
    sum(vals, ret);
    ASSERT((ret == 6), 6, ret);

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_sum2() {
    auto s = NOW;

    mpz_class ret = 0;
    Vec vals = {100,250,333};
    sum(vals, ret);
    ASSERT((ret == 683), 683, ret);

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}


bool test_isbsmooth1() {
    auto s = NOW;
    mpz_class n = 95;
//    mpz_class n = 439;
    mpz_class b = 19;
    IS_TRUE(is_b_smooth(n, b))

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_isbsmooth2() {
    auto s = NOW;
    mpz_class n = 24;
    mpz_class b = 7;
    IS_TRUE(is_b_smooth(n, b))

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_isbsmooth3() {
    auto s = NOW;
    mpz_class n = 24;
    mpz_class b = 2;
    IS_TRUE(!is_b_smooth(n, b))

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_pltb1() {
    auto s = NOW;
    mpz_class b = 30;
    long ret;
    vector<long> m;
    ret = primes_lt_b(b, m);
    ASSERT((ret == 10), 10, ret)

    auto e = NOW;
    auto t = DUR(s, e);

    PASSED;

    return true;
}

bool test_primefacts1() {
    auto s = NOW;
    mpz_class b = 315;
    vector<long> ret = get_factor_vector(b);
    ASSERT((ret[0] == 2), 2, ret[0])
    ASSERT((ret[1] == 1), 1, ret[1])
    ASSERT((ret[2] == 1), 1, ret[2])

    auto e = NOW;
    auto t = DUR(s, e)

    PASSED;

    return true;
}

bool test_evenindices1() {
    auto s = NOW;
    vector<long> v = {2,1,1};
    vector<long> ret = get_even_indices(v, 10);
    ASSERT((ret[0] == 0), 0, ret[0])
    IS_TRUE((ret.size() == 1))

    auto e = NOW;
    auto t = DUR(s, e)

    PASSED;

    return true;
}

bool test_qs1() {
    auto s = NOW;
    mpz_class n, f1, f2;
    n = 315;

    quad_sieve(n, f1, f2);
    cout << f1 << " " << f2 << endl;
    IS_TRUE((f1 * f2 == n))

    auto e = NOW;
    auto t = DUR(s, e)

    PASSED;

    return true;
}


int main() {
    test_prod1();
    test_prod2();
    test_prod3();
    test_sum1();
    test_sum2();
    test_isbsmooth1();
    test_isbsmooth2();
    test_isbsmooth3();
    test_pltb1();
    test_primefacts1();
    test_evenindices1();
    test_qs1();
}
