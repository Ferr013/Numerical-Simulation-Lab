CC = g++
CFLAGS = -Wall -O3 --std=c++11

Annealing.exe : Annealing.o random.o
	$(CC) random.o Annealing.o -o Annealing.exe
Annealing.o : Annealing.cpp
	$(CC) -c Annealing.cpp -o Annealing.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Annealing.exe *.gene evolution.*
