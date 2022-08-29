CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Lez1_11.o function.o random.o
	$(CC) random.o function.o Lez1_11.o -o Lez1_11.exe
Lez1_11.o : Lez1_11.cpp
	$(CC) -c Lez1_11.cpp -o Lez1_11.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
function.o : function.cpp random.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)
clean :
	rm *.o Lez1_11.exe seed.out
