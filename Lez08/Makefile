CC = g++
CFLAGS = -Wall -O3 --std=c++11

Esercizio8.exe : SimAnnealing.o random.o
	$(CC) random.o SimAnnealing.o -o Esercizio8.exe
SimAnnealing.o : SimAnnealing.cpp
	$(CC) -c SimAnnealing.cpp -o SimAnnealing.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Esercizio8.exe
