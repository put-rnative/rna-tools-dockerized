#ifndef Pdb_H
#define Pdb_H

#include "Chain.h"

class Pdb
{
public:
	Pdb();
	Pdb(char *);
	Pdb(string);
	Pdb *copy();
	void readPDB(string);
	void push(Chain *);
	void push(Chain &);
	Chain &operator [](int);
	void updateChains(string);

	void setLen();
	void setResNum();

	int getLen();
	int totalAtoms();
	double getDist(int, int);
	string getSeq();
	
	/* IO function */
	void print();
	void printAsDNA();
	void write(string);
	
	/* assemble function */
	void move(double, double, double);
	void rotate(Matr_ *m);
	void format();
	void addP();
	void mutate(string);
	void rotateByX(double);
	void rotateByZ(double);

	/* attributes */
	string name;
	int len;
	vector<Chain> chains;
};

class Pdbs {
public:
	Pdbs();
	~Pdbs();

	int getLen();
	void resize(int);
	Pdb *at(int);
	Pdb *operator [](int);
	void push(Pdb *);
	
	vector<Pdb *> PdbList;
};

#endif // Pdb_H
