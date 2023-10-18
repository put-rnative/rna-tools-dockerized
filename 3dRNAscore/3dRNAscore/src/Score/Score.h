#ifndef SCORE_H
#define SCORE_H

#include "DistAnal.h"
#include "DihAnal.h"

class Score {
public:
	Score();
	Score(std::string);
	double run(Obj<RNA>);
	void trainDistPar(string);
	void trainDihPar(string);

	string par_dist_obs;
	string par_dist_ref;
	string dihPar;
	string reference_state = "average";
	int cutoff = 20;

	Obj<DistAnal> _distAnal;
	Obj<DihAnal> _dihAnal;

	double _dist_bin = 0.3;
	double _dih_bin = 4.5;

	double _constant = 27.1118;
	double _distWeight = 0.433513;
	double _dihWeight = 1.59348;

	double _distScore = 0;
	double _dihScore = 0;
	double _score = 0;
};



#endif




