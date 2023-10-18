#ifndef DNA_H
#define DNA_H

#include "Chain.h"

class DNA
{
public:
	DNA();
	DNA(char *);
	DNA(string);
	DNA *copy();
	void readPDB(string);
	void push(Chain *);
	void push(Chain &);
	Chain &operator [](int);
	void updateChains(string);

	int setCutOff(int);
	void setLen();

	/* 'get' function */
	int getLen();
	int totalAtoms();
	double getDist(int, int);
	string getSeq();
	
	/* IO function */
	void print();
	void write(string);

	/* attributes */
	string name;
	int len;
	vector<Chain> chains;
};

class DNAs {
public:
	DNAs();
	~DNAs();

	int getLen();
	void resize(int);
	DNA *at(int);
	DNA *operator [](int);
	void push(DNA *);
	
	vector<DNA *> DNAList;
};

#endif // DNA_H_INCLUDED
