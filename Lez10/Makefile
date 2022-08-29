CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

Lez10P.exe : Lez10P.o random.o GenAlgs.o
	$(CC) random.o GenAlgs.o Lez10P.o -o Lez10P.exe
Lez10P.o : Lez10P.cpp GenAlgs.h randomP.h
	$(CC) -c Lez10P.cpp -o Lez10P.o $(CFLAGS)
GenAlgs.o: GenAlgs.h GenAlgs.cpp randomP.h
	$(CC) -c GenAlgs.cpp -o GenAlgs.o $(CFLAGS)
random.o : randomP.cpp randomP.h
	$(CC) -c randomP.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Lez10P.exe
