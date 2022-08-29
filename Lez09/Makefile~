CC = g++
CFLAGS = -Wall -O3 --std=c++11

Lez9.exe : Lez9.o random.o GenAlgs.o
	$(CC) random.o GenAlgs.o Lez9.o -o Lez9.exe
Lez9.o : Lez9.cpp GenAlgs.h random.h
	$(CC) -c Lez9.cpp -o Lez9.o $(CFLAGS)
GenAlgs.o: GenAlgs.h GenAlgs.cpp random.h
	$(CC) -c GenAlgs.cpp -o GenAlgs.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Lez9.exe
