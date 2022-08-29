CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Lez3.o function.o random.o Funzioni.o
	$(CC) random.o function.o  Funzioni.o Lez3.o -o Lez3.exe
Lez3.o : Lez3.cpp
	$(CC) -c Lez3.cpp -o Lez3.o $(CFLAGS)
random.o : random.cpp random.h Funzioni.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
function.o : function.cpp random.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)
Funzioni.o : Funzioni.cpp random.h
	$(CC) -c Funzioni.cpp -o Funzioni.o $(CFLAGS)
clean :
	rm *.o Lez3.exe seed.out
