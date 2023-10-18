CC = g++

all: matrix.o string.o
clean: clear

matrix.o: matrix.cpp matrix.h
	${CC} ${CFLAG} -c matrix.cpp -o matrix.o

string.o: string.cpp string.h
	${CC} ${CFLAG} -c string.cpp -o string.o

clear:
	\rm -f *.o
