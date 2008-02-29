all:
	gcc -pedantic -Wall -O3 -march=nocona -c *.c
	gcc -lgmp -O3 -march=nocona -o tester *.o

test:
	gcc -g -Wall -pedantic -std=c99 -c scalar.c pol.c test.c
	gcc -g -o tester scalar.o pol.o test.o

clean:
	rm -f tester tijdelijk gmon.out
	rm -f *.o

debug:
	gcc -g -DKIJKEN -Wall -c *.c
	gcc -lgmp -g -Wall -o tester *.o

profiler:
	gcc -pg -DPROFILER -Wall -O2 -march=nocona -c *.c
	gcc -pg -Wall -lgmp -O2 -march=nocona -o tester *.o
