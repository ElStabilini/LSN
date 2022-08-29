CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Lez2_2b.o function.o random.o
	$(CC) random.o function.o Lez2_2b.o -o Lez2_2b.exe
Lez2_2b.o : Lez2_2b.cpp
	$(CC) -c Lez2_2b.cpp -o Lez2_2b.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
function.o : function.cpp random.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)
RandomWalk.o : RandomWalk.cpp random.h
	$(CC) -c RandomWalk.cpp -o RandomWalk.o $(CFLAGS)
clean :
	rm *.o Lez2_2b.exe seed.out
