#include <gmpxx.h>

using namespace std;


bool is_b_smooth(mpz_class n, mpz_class &b);
void quad_sieve(const mpz_class &n, mpz_class &fact1, mpz_class &fact2);
long long primes_lt_b(mpz_class &b);