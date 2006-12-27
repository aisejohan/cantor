all:
	gcc -pedantic -Wall -O3 -march=nocona -c *.c
	gcc -lgmp -O3 -march=nocona -o tester *.o

clean:
	rm -f tester tijdelijk gmon.out
	rm -f *.o

debug:
	gcc -g -DKIJKEN -Wall -c *.c
	gcc -lgmp -g -Wall -o tester *.o

profiler:
	gcc -g -pg -DPROFILER -Wall -c *.c
	gcc -g -pg -Wall -lgmp -o tester *.o
