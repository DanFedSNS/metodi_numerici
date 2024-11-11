CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99 

binaries=mcmc_analysis0 mcmc_analysis1 mcmc_analysis2 mcmc_analysis3 ising2d_metro \
	ising_metro ising_analysis ising_wolff on_metro on_wolff

all: $(binaries)

random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
	$(CC) $(CFLAGS) -c ../lib/random.c

pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
	$(CC) $(CFLAGS) -c ../lib/pcg32min.c

boxmuller.o: ../lib/boxmuller.c ../include/boxmuller.h
	$(CC) $(CFLAGS) -c ../lib/boxmuller.c

read_data.o: ../lib/read_data.c ../include/read_data.h
	$(CC) $(CFLAGS) -c ../lib/read_data.c

mcmc_analysis0: ../src/mcmc_analysis0.c read_data.o
	$(CC) $(CFLAGS) -c ../src/mcmc_analysis0.c 
	$(CC) $(CFLAGS) mcmc_analysis0.o read_data.o -o $@ -lm

mcmc_analysis1: ../src/mcmc_analysis1.c read_data.o
	$(CC) $(CFLAGS) -c ../src/mcmc_analysis1.c 
	$(CC) $(CFLAGS) mcmc_analysis1.o read_data.o -o $@ -lm

mcmc_analysis2: ../src/mcmc_analysis2.c pcg32min.o random.o read_data.o
	$(CC) $(CFLAGS) -c ../src/mcmc_analysis2.c  
	$(CC) $(CFLAGS) mcmc_analysis2.o pcg32min.o random.o read_data.o -o $@ -lm

mcmc_analysis3: ../src/mcmc_analysis3.c read_data.o
	$(CC) $(CFLAGS) -c ../src/mcmc_analysis3.c  
	$(CC) $(CFLAGS) mcmc_analysis3.o read_data.o -o $@ -lm

geometry.o: ../lib/geometry.c ../include/geometry.h 
	$(CC) $(CFLAGS) -c ../lib/geometry.c  

ising2d_metro: ../src/ising2d_metro.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/ising2d_metro.c  
	$(CC) $(CFLAGS) ising2d_metro.o random.o pcg32min.o -o $@ -lm

ising3d_metro: ../src/ising3d_metro.c random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/ising3d_metro.c  
	$(CC) $(CFLAGS) ising3d_metro.o random.o pcg32min.o -o $@ -lm

ising_metro: ../src/ising_metro.c geometry.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/ising_metro.c  
	$(CC) $(CFLAGS) geometry.o ising_metro.o random.o pcg32min.o -o $@ -lm

ising_analysis: ../src/ising_analysis.c read_data.o
	$(CC) $(CFLAGS) -c ../src/ising_analysis.c  
	$(CC) $(CFLAGS) ising_analysis.o read_data.o -o $@ -lm

ising_wolff: ../src/ising_wolff.c geometry.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/ising_wolff.c  
	$(CC) $(CFLAGS) geometry.o ising_wolff.o random.o pcg32min.o -o $@ -lm

nvector.o: ../lib/nvector.c ../include/nvector.h ../include/random.h 
	$(CC) $(CFLAGS) -c ../lib/nvector.c  

on_metro: ../src/on_metro.c geometry.o nvector.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/on_metro.c  
	$(CC) $(CFLAGS) geometry.o nvector.o on_metro.o random.o pcg32min.o -o $@ -lm

on_wolff: ../src/on_wolff.c geometry.o nvector.o random.o pcg32min.o
	$(CC) $(CFLAGS) -c ../src/on_wolff.c  
	$(CC) $(CFLAGS) geometry.o nvector.o on_wolff.o random.o pcg32min.o -o $@ -lm

.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
