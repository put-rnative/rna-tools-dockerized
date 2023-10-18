#ifndef RESIDUE_H_INCLUDED
#define RESIDUE_H_INCLUDED

#include "Atom.h"

class Residue
{
public:
	Residue();
	Residue(vector<string> &, string);
	Residue(Point *, int, int);
	Atom &operator [](int);
	Point *getBaseMassCenter();
	int getAtomPos(string, Point &);
	int getAmount();
	Point *getBaseVec();
	int nextTo(Residue &);

	int num;
	string number;
	int atomNum;
	string rnaName;
	string name;
	vector<Atom> atoms;
};

#endif // RESIDUE_H_INCLUDED
