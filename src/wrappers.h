#pragma once
#include <gmpxx.h>

inline void powm(mpz_class &res, const mpz_class &base, const mpz_class &exp, const mpz_class &mod) {
    mpz_powm(res.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), mod.get_mpz_t());
}
inline void powm_ui(mpz_class &res, const mpz_class &base, const ulong exp, const mpz_class &mod) {
    mpz_powm_ui(res.get_mpz_t(), base.get_mpz_t(), exp, mod.get_mpz_t());
} 

inline void pow_ui(mpz_class &res, const mpz_class &base, const ulong exp) {
    mpz_pow_ui(res.get_mpz_t(), base.get_mpz_t(), exp);
}

inline int tstbit(const mpz_class &a, ulong bit) {
    return mpz_tstbit(a.get_mpz_t(), bit);
}

inline int legendre(const mpz_class &n, const mpz_class &p) {
    mpz_class res;
    powm(res, n, (p- 1) / 2, p);
    return res.get_si();
}

inline void invert(mpz_class &res, const mpz_class &a, const mpz_class &b) {
    mpz_invert(res.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}

inline void fdiv_q_2exp(mpz_class &res, const mpz_class &a, ulong exp) {
    mpz_fdiv_q_2exp(res.get_mpz_t(), a.get_mpz_t(), exp);
}

inline void sqrt(mpz_class &n) {
    mpz_sqrt(n.get_mpz_t(), n.get_mpz_t());
}

inline mpz_class safe_mod(mpz_class a, mpz_class b) {
    return (a % b + b) % b;
}

inline int digit_length(const mpz_class &n, int base) {
    return mpz_sizeinbase(n.get_mpz_t(), base);
}