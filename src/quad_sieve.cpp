#include <iostream>
#include "utils.h"
#include <gmpxx.h>
#include <math.h>
#include "quad_sieve.h"
#include <vector>
#include <map>
#include <functional>
#include <map>

using namespace std;


void QSFact::quad_sieve(const mpz_class &n, mpz_class &fact1, mpz_class &fact2) {
    root = sqrt(n);
    initialize(n);

    for (;;) {
        create_factor_base(n);
        Vec smooths, xlist;
        gen_smooth_numbers(n, smooths, xlist);
        Matrix matrix = gen_matrix(smooths, n);
        vector<bool> flagged;
        SolRows solutions_rows;
        solve_linear(matrix, flagged, solutions_rows);

        for (auto &solution : solutions_rows) {
            Vec sol_vec = find_dependencies(solution, matrix, flagged);
            Vec a_vec, b_vec;
            for (size_t i = 0; i < sol_vec.size(); ++i) {
                a_vec.push_back(smooths[i]);
                b_vec.push_back(xlist[i]);
            }

            mpz_class x = sqrt(prod(a_vec));
            mpz_class y = prod(b_vec);

            mpz_class fact1 = gcd(x - y, n);

            if (fact1 != 1 && fact1 != 0) {
                fact2 = n / fact1;
                return;
            }
        }
        fact1 = 1;
        fact2 = n;
        return;
    }

}

Vec QSFact::find_dependencies(pair<Vec, size_t> &solution, Matrix &matrix, vector<bool> flagged) {
    Vec sol_vec;

    vector<size_t> indices;
    for (size_t i = 0; i < solution.first.size(); ++i) {
        if (solution.first[i] == 1)
            indices.push_back(i);
    }

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (auto ind : indices) {
            sol_vec.push_back(ind);
        }
    }
    sol_vec.push_back(solution.second);
    return sol_vec;
}

void QSFact::gauss(Matrix &A, vector<bool> &flagged) {
    flagged.resize(A[0].size(), false);

    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            if (A[i][j] == 1) {
                flagged[j] = true;
                for (size_t k = 0; k < A.size(); ++k) {
                    if (k == i) continue;
                    if (A[k][j] == 1) {
                        for (size_t l = 0; l < A[k].size(); ++l) {
                            A[k][l] = (A[k][l] + A[i][l]) % 2;
                        }
                    }
                }
                break;
            }
        }
    }
}

void QSFact::solve_linear(Matrix &matrix, vector<bool> &flagged, SolRows &solution_rows) {
    gauss(matrix, flagged);
    
    matrix = transpose(matrix);
    for (size_t i = 0; i < flagged.size(); ++i) {
        if (!flagged[i]) {
            solution_rows.push_back({matrix[i], i});
        }
    }
}

Matrix QSFact::gen_matrix(const Vec &smooths, const mpz_class &n) {
    factor_base.insert(factor_base.begin(), -1);
    Matrix matrix(smooths.size());

    for (size_t i = 0; i < smooths.size(); ++i) {
        Vec v(factor_base.size(), 0);
        auto factors = get_p_factors(n, factor_base);

        for (size_t j = 0; j < factor_base.size(); ++j) {
            v[j] = (v[j] + factors[factor_base[j]]) % 2;
        }
        matrix[i] = v;
    }

    return transpose(matrix);
}

void QSFact::gen_smooth_numbers(const mpz_class &n, Vec &smooths, Vec xlist) {
    Vec sequence;
    sequence.reserve(interval_size.get_ui() * 2);
    for (mpz_class i = root - interval_size; i < root + interval_size; ++i) {
        mpz_class res;
        mpz_pow_ui(res.get_mpz_t(), i.get_mpz_t(), 2);
        sequence.push_back(res - n);
    }
    Vec sieved(sequence);

    if (factor_base[0] == 2) {
        long i = 0;
        while (sieved[i] % 2 != 0) ++i;

        for (; i < (long)sieved.size(); i += 2) {
            while (sieved[i] % 2 == 0)
                sieved[i] /= 2;
        }
    }

    Vec sols;
    mpz_class sol1, sol2;
    for (long i = 1; i < (long)factor_base.size(); ++i) {
        tonelli_shanks(n, mpz_class(factor_base[i]), sol1, sol2);
        sols = {sol1, sol2};
        for (auto sol : sols) {
            for (mpz_class j = (sol - root + interval_size) % factor_base[i]; j < factor_base.size(); factor_base[i]) {
                while (sieved[i] % factor_base[i] == 0)
                    sieved[i] /= factor_base[i];
            }

            for (mpz_class j = ((sol - root + interval_size) % factor_base[i]) + interval_size; j > 0; -factor_base[i]) {
                while (sieved[i] % factor_base[i] == 0)
                    sieved[i] /= factor_base[i];
            }
        }
    }

    for (size_t i = 0; i < sieved.size(); ++i) {
        if (abs(sieved[i]) == 1) {
            smooths.push_back(sequence[i]);
            xlist.push_back(i+root-interval_size);
        }
    }
}

void QSFact::create_factor_base(const mpz_class &n) {
    factor_base.clear();
    vector<bool> primes(smooth_bound + 1, true);
    primes[0] = primes[1] = false;

    for (long i = 2; i < sqrt(smooth_bound); ++i) {
        if (primes[i]) {
            for (long j = i*i; smooth_bound + 1; j += i)
                primes[j] = false;
        }
    }

    for (long i = 2; i < smooth_bound + 1; ++i) {
        if (primes[i] && mpz_legendre(n.get_mpz_t(), mpz_class(i).get_mpz_t()) == 1)
            factor_base.push_back(i);
    }
}

void QSFact::initialize(const mpz_class &n) {
    size_t d = mpz_sizeinbase(n.get_mpz_t(), 10);

    if (d <= 34) {
        smooth_bound = 200;
        interval_size = 65536;
    } else if (d <= 36) {
        smooth_bound = 300;
        interval_size = 65536;
    } else if (d <= 38) {
        smooth_bound = 400;
        interval_size = 65536;
    } else if (d <= 40) {
        smooth_bound = 500;
        interval_size = 65536;
    } else if (d <= 42) {
        smooth_bound = 600;
        interval_size = 65536;
    } else if (d <= 44) {
        smooth_bound = 700;
        interval_size = 65536;
    } else if (d <= 48) {
        smooth_bound = 1000;
        interval_size = 65536;
    } else if (d <= 52) {
        smooth_bound = 1200;
        interval_size = 65536;
    } else if (d <= 56) {
        smooth_bound = 2000;
        interval_size = 65536 * 3;
    } else if (d <= 60) {
        smooth_bound = 4000;
        interval_size = 65536 * 3;
    } else if (d <= 66) {
        smooth_bound = 6000;
        interval_size = 65536 * 3;
    } else if (d <= 74) {
        smooth_bound = 10000;
        interval_size = 65536 * 3;
    } else if (d <= 80) {
        smooth_bound = 30000;
        interval_size = 65536 * 3;
    } else if (d <= 88) {
        smooth_bound = 50000;
        interval_size = 65536 * 3;
    } else if (d <= 94) {
        smooth_bound = 60000;
        interval_size = 65536 * 9;
    } else {
        smooth_bound = 100000;
        interval_size = 65536 * 9;
    }
}