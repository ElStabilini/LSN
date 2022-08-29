CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Lez1_3.o function.o random.o
	$(CC) random.o function.o Lez1_3.o -o Lez1_3.exe
Lez1_3.o : Lez1_3.cpp
	$(CC) -c Lez1_3.cpp -o Lez1_3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
function.o : function.cpp random.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)
clean :
	rm *.o Lez1_3.exe seed.out
