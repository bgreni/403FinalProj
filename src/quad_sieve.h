#include <gmpxx.h>

using namespace std;

void quad_sieve(const mpz_class &n, mpz_class &fact1, mpz_class &fact2);
long primes_lt_b(mpz_class &b, vector<long> &prime_set);
vector<long> get_factor_vector(const mpz_class n);