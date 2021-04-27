#include <gmpxx.h>
#include <iostream>
#include "quad_sieve.h"
#include "utils.h"
#include <assert.h>
#include "test_numbers.h"
#include <map>

using namespace std;

map<int, mpz_class> difficulties {
    {1, RSA_32bit},
    {2, RSA_45bit},
    {3, RSA_60bit},
    {4, RSA_64bit},
    {5, RSA_80bit},
    {6, RSA_80bit_2},
    {7, RSA_90bit},
    {8, RSA_100bit},
    {9, RSA_129bit},
};


int main(int argc, char* argv[]) {

    if (argc != 3) {
        cout << "requires 2 args, recieved: " << argc-1 << endl;
        return 0;
    }

    string arg1(argv[1]);
    string arg2(argv[2]);
    int level = 0;
    mpz_class n;
    if (arg1 == "--level") {
        level = stoi(arg2);
        if (level < 0 || level > 9)  {
            cout << RED << "Error: level must with the range [0, 9], recieved " << level << RESET << endl;
            return 0;
        }
        n = difficulties[level];
    } else if (arg1 == "--user") {
        n = arg2;
    } else {
        cout << RED << "Error: invalid arg given\nfirst arg must be --level or --user\n" << arg2 << " recieved" << RESET << endl;
        return 0;
    }
    

    if (level > 6 || digit_length(n, 2) > 80)  {
        char resp;
        cout << YELLOW << "Warning: factoring integer over roughly 80 bits may take a very long time" << RESET << endl;
        cout << "Do you still wish to continue? [Y/n] ";
        cin >> resp;
        if (tolower(resp) != 'y')
            return 1;
    }

    cout << "Attempting to factor " << n << "..." << endl;
    mpz_class f1, f2;;
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