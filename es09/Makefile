CC = g++
CFLAGS = -Wall -O3 --std=c++11

GA_TSP.exe : GA_TSP.o random.o
	$(CC) random.o GA_TSP.o -o GA_TSP.exe
GA_TSP.o : GA_TSP.cpp
	$(CC) -c GA_TSP.cpp -o GA_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o GA_TSP.exe *.gene evolution.dat
