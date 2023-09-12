#ifndef RNA_H
#define RNA_H

#include "Chain.h"

class RNA
{
public:
	RNA();
	RNA(char *);
	RNA(string);
	RNA *copy();
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

class RNAs {
public:
	RNAs();
	~RNAs();

	int getLen();
	void resize(int);
	RNA *at(int);
	RNA *operator [](int);
	void push(RNA *);
	
	vector<RNA *> RNAList;
};

#endif // RNA_H
