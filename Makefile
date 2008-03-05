all:
	gcc -pedantic -Wall -O3 -march=nocona -c *.c
	gcc -O3 -march=nocona -o tester *.o

test:
	gcc -g -Wall -pedantic -std=c99 -c scalar.c pol.c test.c
	gcc -g -o tester scalar.o pol.o test.o

clean:
	rm -f tester tijdelijk gmon.out
	rm -f *.o

debug:
	gcc -g -DKIJKEN -Wall -c *.c
	gcc -g -Wall -o tester *.o

profiler:
	gcc -pg -DPROFILER -Wall -O2 -march=nocona -c *.c
	gcc -pg -Wall -O2 -march=nocona -o tester *.o
