#ifndef CHAIN_H_INCLUDED
#define CHAIN_H_INCLUDED

#include "Residue.h"

class Chain
{
public:
	Chain();
	Chain(vector<string> &, string);

	void push(Residue *);
	void push(Residue &);
	Residue &operator [](int);

	string name;
	string rnaName;
	vector<Residue> residues;
};

#endif // CHAIN_H_INCLUDED
