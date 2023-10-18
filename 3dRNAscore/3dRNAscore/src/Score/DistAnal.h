#ifndef DISTANAL_H
#define DISTANAL_H

#include "../Pdb.h"

class DistAnal {
public:
	DistAnal(int = 20, double = 0.5, string = "average");
	~DistAnal();

//	void readRNA(RNA *);
	void readRNA(Obj<RNA>);
	void train();
	double scoring();

	void readObsParm(string = "");
	void readRefParm(string = "");
	void initObsProb();
	void initRefProb();

	void printObsParm();
	void printRefParm();
	void printObsProb();
	void printRefProb();

	double getScore();

private:
	string reference_state;
	int *num;
	int *type;
	int *ntLen;
	Point **list;
	int len;

	int *obsParm;
	int *refParm;
	double *obsProb;
	double *refProb;

	double score;

	double interval;
	int cutoff;
	int bins;
	double penalty;
};

#endif

