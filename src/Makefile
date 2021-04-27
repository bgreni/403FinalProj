COMP = g++ -Wall -std=c++17 -Ofast -g -o
OBJ = g++ -Wall -std=c++17 -Ofast -g -c


defaults: main

main: main.o quad_sieve.o utils.o
	$(COMP) main main.o quad_sieve.o utils.o -lgmpxx -lgmp

tests: tests.o quad_sieve.o utils.o
	$(COMP) tests tests.o quad_sieve.o utils.o -lgmpxx -lgmp

%.o: %.cpp
	$(OBJ) $< -lgmpxx -lgmp

main.o: main.cpp utils.h quad_sieve.h test_numbers.h
tests.o: tests.cpp utils.h quad_sieve.h wrappers.h test_numbers.h
quad_sieve.o: quad_sieve.cpp quad_sieve.h utils.h wrappers.h
utils.o: utils.cpp utils.h wrappers.h

clean:
	-rm *.o $(objects)
	-rm main
	-rm tests
	-rm vgcore*