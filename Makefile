all:
	export OMP_DYNAMIC=true
	gcc -o computeOOP computeOOP.c -fopenmp -lm -Wall -O3
	./computeOOP
