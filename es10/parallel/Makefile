CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

MPI_Annealing.exe : MPI_Annealing.o random.o
	$(CC) random.o MPI_Annealing.o -o MPI_Annealing.exe
MPI_Annealing.o : MPI_Annealing.cpp
	$(CC) -c MPI_Annealing.cpp -o MPI_Annealing.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o MPI_Annealing.exe *.gene *evo*