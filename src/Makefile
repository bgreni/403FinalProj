LINKS = -lgmpxx -lgmp
COMP = g++ $(LINKS) -Wall -std=c++17 -O3 -o
OBJ = g++ $(LINKS) -Wall -std=c++17 -O3 -c


defaults: qsmain

qsmain: qsmain.o quad_sieve.o utils.o
	$(COMP) qsmain qsmain.o quad_sieve.o utils.o

tests: tests.o quad_sieve.o utils.o
	$(COMP) tests tests.o quad_sieve.o utils.o

%.o: %.cpp
	$(OBJ) $<

qsmain.o: qsmain.cpp utils.h quad_sieve.h test_numbers.h
tests.o: tests.cpp utils.h quad_sieve.h wrappers.h test_numbers.h pollard.h
quad_sieve.o: quad_sieve.cpp quad_sieve.h utils.h wrappers.h
utils.o: utils.cpp utils.h wrappers.h

clean:
	-rm *.o $(objects)
	-rm main
	-rm tests
	-rm vgcore*