#ifndef ATOM_H_INCLUDED
#define ATOM_H_INCLUDED

#include "../Tools.h"

class Atom {
public:
	Atom();
	Atom(string &, string);
	Atom(Point, string, string, int);
	Point *coord();
	double dist(Atom &);

	string name;
	string resName;
	string rnaName;
	string line;
	int num;
	double mass;
	double x, y, z;
	int flag;
};

#endif // ATOM_H_INCLUDED
