all:
	gcc -pedantic -Wall -O3 -march=nocona -c main.c utils.c pol.c scalar.c xu_and_sparse.c list_degrees.c
	gcc -O3 -march=nocona -o tester main.o utils.o pol.o scalar.o xu_and_sparse.o list_degrees.c

test:
	gcc -g -Wall -pedantic -std=c99 -c utils.c list_degrees.c pol.c scalar.c test.c  xu_and_sparse.c
	gcc -g -o tester utils.o list_degrees.o pol.o scalar.o test.o xu_and_sparse.o

clean:
	rm -f tester tijdelijk gmon.out
	rm -f *.o

debug:
	gcc -g -DKIJKEN -Wall -c main.c pol.c scalar.c xu_and_sparse.c
	gcc -g -Wall -o tester main.o pol.o scalar.o xu_and_sparse.o

profiler:
	gcc -pg -O2 -Wall -march=nocona -c main.c utils.c list_degrees.c pol.c scalar.c xu_and_sparse.c
	gcc -pg -O2 -Wall -march=nocona -o tester utils.c list_degrees.c main.o pol.o scalar.o xu_and_sparse.o

gcov:
	gcc -fprofile-arcs -ftest-coverage -Wall -march=nocona -c main.c pol.c scalar.c xu_and_sparse.c
	gcc -fprofile-arcs -ftest-coverage -march=nocona -o tester main.o pol.o scalar.o xu_and_sparse.o
