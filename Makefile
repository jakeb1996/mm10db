CC = gcc
NEW_CC = gcc-4.9

FLAGS = -Wall -pedantic -ansi -O3
LIB_OPENMP = -fopenmp 
LIB_THREAD = -pthread

EXE_OPENMP = findMismatches_openMP
EXE_THREAD = findMismatches_threads

all: openMP pthread

openMP: findMismatches_openMP.c 
	$(NEW_CC) $(LIB_OPENMP) $(FLAGS) findMismatches_openMP.c -o $(EXE_OPENMP)

pthread: findMismatches_threads.c
	$(CC) $(LIB_THREAD) $(FLAGS) findMismatches_threads.c -o $(EXE_THREAD)

clean:
	rm $(EXE_THREAD) $(EXE_OPENMP)
