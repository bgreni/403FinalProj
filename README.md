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
1. Enter the source folder: `cd src`

### Running the main program
2. compile the project: `make`
3. Run the program with the appropriate arguments

    For choose the number to factor, there are two options, use `--level` and choose a
    number [1, 9]. Each level corresponds to a RSA number generated [here](https://bigprimes.org/RSA-challenge), with difficulty of the number increasing with each level. The exact numbers the levels refer to can be found in `qsmain.cpp`
    ```
    ./qsmain --level 5
    ```
    Another option is to input a number yourself using `--user`
    ```
    ./qsmain --user 23894623598
    ```
    An additional optional argument is `--iter-cap` which allows the user to set the
    max number of iterations the algorithm will try before giving up (default: 100).
    This should be used if you want to test out a large number, but don't want to
    risk it running for hours.
    ```
    ./qsmain --level 8 --iter-cap 20
    ```
### Running Tests
4. Running Unit tests
    
    A few unit tests have also been included that includes test for number up until 80 bits long, as well as one test which tests random number with factors from with the interval [400, sqrt(MAX_INT)]. The random test will also report the number of "weak factors" that were produced. I have defined a weak factor as anything less than 10. This metric isn't terribly important, but since the factor of every test case is at least 400, then it should have been possible to find a more interesting factor.
    ```
    make tests
    ./tests
    ```

## Program Ouput
The output is simply the time it took to find a solution in seconds, as well as
the two factors that were found, e.g:
```
Attempting to factor 2345346...
time taken (seconds): 0.0891863
factor 1: 8986
factor 2: 261
```

## Issues to Note
Very large inputs (~90 bits) seem to have a chance of causing a segfault, and I was never able to understand why, as running the program through valgrind is orders of magnitude slower, and it may take far longer to find the segfault than I am at liberty to wait. It is also possible that the random tests section of the unit tests may take far longer than the average of 4 seconds, if a problematic number is generated, of which I also don't know the cause. Otherwise, the program should work as intended.

Another note is this project was developed using my personal machine with these specs:
```
Ryzen 3700x 8-core CPU
32 GB RAM
Manjaro Linux
gcc version 10.2.0
```
I seem to have ironed out any MacOS-specific issues due to the differences in the C compilers on Linux and MacOS, but just in case something goes wrong, I would first try and run it on a linux system.

## Summary of Files
- qsmain.cpp: The entry point for the CLI program, usage explained [above](##Running-the-program)
- quad_sieve.cpp: Contains the core implementation of the SIQS algorithm
- utils.cpp/.h: A collection of utilities for debugging, more general algorithms e.g: matrix transpose, and helpful macros and type definitions
- pollard.h: An implementation of Pollard Rho's to compare the performance of my SIQS algorithm
- wrappers.h: A collection of wrapper functions for the old C-style GMP functions, so I could avoid having to use `mpz_class.get_mpz_t()` so often
- test_numbers.h: All the numbers I used for unit testing, as well as a few that are too large for me to factor in a reasonable timeframe

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

 self intializing values taken from https://github.com/skollmann/PyFactorise/blob/master/factorise.py

all RSA numbers generated here https://bigprimes.org/RSA-challenge

https://stackoverflow.com/questions/9158150/colored-output-in-c/9158263