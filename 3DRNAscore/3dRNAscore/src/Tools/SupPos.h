#ifndef SUPPOS_H
#define SUPPOS_H

#include "std.h"
#include "Point.h"
#include "Matr_.h"

class SupPos {
public:
	SupPos(Matr_ *, Matr_ *);
	~SupPos();
	Matr_ *getRot();
	Point *getC1();
	Point *getC2();
	double getRMSD();

private:
	int len;
	Matr_ *rot;
	Point c1;
	Point c2;
	double rmsd;
};



#endif






