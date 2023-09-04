CC = g++

TOOLIB = ../tool-library

all: main
clean: clear

main: pdb.o model.o chain.o residue.o atom.o format.o

pdb.o: pdb.cpp pdb.h format.h ${TOOLIB}/string.h
	${CC} ${CFLAG} -c pdb.cpp -o pdb.o 

model.o: model.cpp model.h format.h
	${CC} ${CFLAG} -c model.cpp -o model.o 

chain.o: chain.cpp chain.h format.h
	${CC} ${CFLAG} -c chain.cpp -o chain.o 

atom.o: atom.cpp atom.h point.h index.h ${TOOLIB}/string.h
	${CC} ${CFLAG} -c atom.cpp -o atom.o 

residue.o: residue.cpp residue.h atom.h point.h index.h  format.h ${TOOLIB}/string.h
	${CC} ${CFLAG} -c residue.cpp -o residue.o 

format.o: format.cpp format.h
	${CC} ${CFLAG} -c format.cpp -o format.o 

clear:
	\rm -f *.o
