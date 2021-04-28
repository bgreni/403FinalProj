# QSFact

## The algorithm

For my project I implemented a self-initializing quadratic sieve
(SIQS). The purpose of this algorithm is to factor extremely large
integers (greater than 10^50). It is the 2 fastest classical integer factoring algorithm, behind the general number field sieve, but SIQS is much simpler.

### Runtime
according to [this](https://mathworld.wolfram.com/QuadraticSieve.html) page from wolfram mathworld, general runtime of the algorithm is O(X) where X = exp(sqrt(ln(n) ln(ln(n)))). From what I
understand, this is due to the fact that the analytically optimal choice of smoothness
B is X ^ (sqr(2)/4) [reference (page 4)](http://www.damianball.com/pdf/portfolio/quadratic-sieve.pdf), and the asymptotic runtime of the algorithm is relative to the size of the
factor base. My program does not use that formula however, as I found computing
logs using gmp is not supported, and I didn't like any of the workarounds I found.
So instead I use a set of fixed parameters based on the base 10 digit length of N. Therefore, my solution is likely some factor slower than optimal.

## Running the program
1. Enter the source folder `cd`


# References
https://en.wikipedia.org/wiki/Quadratic_sieve

https://en.wikipedia.org/wiki/Smooth_number#Definition

https://planetmath.org/QuadraticSieve

http://www.damianball.com/pdf/portfolio/quadratic-sieve.pdf

https://www.geeksforgeeks.org/p-smooth-numbers-p-friable-number/

https://www.youtube.com/watch?v=Y3N0vZoPCWE

https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

https://martinlauridsen.info/pub/bsc_thesis.pdf

https://github.com/Maosef/Quadratic-Sieve/blob/242238b74ded1983f91160952e57e949beae16eb/Quadratic%20Sieve.py

https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm

https://github.com/cheran-senthil/PyRival/blob/master/pyrival/algebra/mod_sqrt.py

self intializing values taken from

https://github.com/skollmann/PyFactorise/blob/master/factorise.py

all RSA numbers generated here

https://bigprimes.org/RSA-challenge

https://stackoverflow.com/questions/9158150/colored-output-in-c/9158263