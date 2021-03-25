COMP = g++ -Wall -std=c++17 -Ofast -o
OBJ = g++ -Wall -std=c++17 -Ofast -c


defaults: main

main: main.o quad_sieve.o
	$(COMP) main main.o quad_sieve.o -lgmpxx -lgmp

tests: tests.o quad_sieve.o
	$(COMP) tests tests.o quad_sieve.o -lgmpxx -lgmp

%.o: %.cpp
	$(OBJ) $< -lgmpxx -lgmp

main.o: main.cpp utils.h quad_sieve.h
tests.o: tests.cpp utils.h quad_sieve.h
quad_sieve.o: quad_sieve.cpp quad_sieve.h utils.h

clean:
	-rm *.o $(objects)
	-rm main
	-rm tests