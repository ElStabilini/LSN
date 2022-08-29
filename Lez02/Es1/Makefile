CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Lez2_1.o function.o random.o Funzioni.o Integrazione.o
	$(CC) random.o function.o  Funzioni.o Integrazione.o Lez2_1.o -o Lez2_1.exe
Lez2_1.o : Lez2_1.cpp
	$(CC) -c Lez2_1.cpp -o Lez2_1.o $(CFLAGS)
random.o : random.cpp random.h Funzioni.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
function.o : function.cpp random.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)
Funzioni.o : Funzioni.cpp random.h
	$(CC) -c Funzioni.cpp -o Funzioni.o $(CFLAGS)
Integrazione.o : Integrazione.cpp Integrazione.h random.h Funzioni.h
	$(CC) -c Integrazione.cpp -o Integrazione.o $(CFLAGS)
clean :
	rm *.o Lez2_1.exe seed.out
