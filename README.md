# QSFact

## The algorithm

For my project I implemented a self-initializing quadratic sieve
(SIQS). The purpose of this algorithm is to factor extremely large
integers (greater than 10^50). It is the 2 fastest classical integer factoring algorithm, behind the general number field sieve, but SIQS is much simpler.

### Runtime
according to [this](https://mathworld.wolfram.com/QuadraticSieve.html) page from wolfram mathworld, general runtime of the algorithm is O(exp(sqrt(ln(n) ln(ln(n))))), and frankly, I have absolutely no idea why. However, I assume that bound is for one iteration of the algorithm, and since my implementation allows for retry iterations if the first attemp fails, then the a more precise expression would be O(exp(sqrt(ln(n) ln(ln(n)))) * NUM_ITERATIONS)

## Running the program



# References
https://en.wikipedia.org/wiki/Quadratic_sieve

https://en.wikipedia.org/wiki/Smooth_number#Definition

https://planetmath.org/QuadraticSieve

http://www.damianball.com/pdf/portfolio/quadratic-sieve.pdf

https://www.geeksforgeeks.org/p-smooth-numbers-p-friable-number/

https://www.youtube.com/watch?v=Y3N0vZoPCWE

https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

https://primes.utm.edu/howmany.html

https://people.cs.clemson.edu/~goddard/MINI/2004/BowmanCochran.pdf

https://martinlauridsen.info/pub/bsc_thesis.pdf

https://github.com/Maosef/Quadratic-Sieve/blob/242238b74ded1983f91160952e57e949beae16eb/Quadratic%20Sieve.py

https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm

https://github.com/cheran-senthil/PyRival/blob/master/pyrival/algebra/mod_sqrt.py

self intializing values taken from

https://github.com/skollmann/PyFactorise/blob/master/factorise.py

all RSA numbers generated here

https://bigprimes.org/RSA-challenge

https://stackoverflow.com/questions/9158150/colored-output-in-c/9158263