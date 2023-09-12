#include "Score.h"

Score::Score() {
	string lib = getenv("RNAscore");
	par_dist_obs = lib + "/lib/par_dist";
	par_dist_ref = "";
	dihPar = lib + "/lib/par_dih";

	_distAnal = new DistAnal(cutoff, _dist_bin, reference_state);
	_dihAnal = new DihAnal(_dih_bin);
	_distAnal->readObsParm(par_dist_obs);
	_distAnal->readRefParm(par_dist_ref);
	_dihAnal->readParm(dihPar);
}

Score::Score(std::string parm) {
	string lib = getenv("RNAscore");
	par_dist_obs = lib + "/lib/par_dist";
	par_dist_ref = "";
	dihPar = lib + "/lib/par_dih";

	ifstream ifile(parm.c_str());
	while (ifile) {
		string line;
		getline(ifile, line);
		vector<string> splited_line;
		tokenize(line, splited_line, " ,:");
		if (splited_line.size() != 2) continue;
		if (splited_line[0] == "par_dist") {
			par_dist_obs = splited_line[1];
		} else if (splited_line[0] == "par_dih") {
			dihPar = splited_line[1]; 
		} else if (splited_line[0] == "cutoff") {
			cutoff = atoi(splited_line[1].c_str());
		} else if (splited_line[0] == "dist_bin_width") {
			_dist_bin = atof(splited_line[1].c_str());
		} else if (splited_line[0] == "dih_bin_width") {
			_dih_bin = atof(splited_line[1].c_str());
		} else if (splited_line[0] == "dist_weight") {
			_distWeight = atof(splited_line[1].c_str());
		} else if (splited_line[0] == "dih_weight") {
			_dihWeight = atof(splited_line[1].c_str());
		}
	}
	ifile.close();

	_distAnal = new DistAnal(cutoff, _dist_bin, reference_state);
	_dihAnal = new DihAnal(_dih_bin);
	_distAnal->readObsParm(par_dist_obs);
	_distAnal->readRefParm(par_dist_ref);
	_dihAnal->readParm(dihPar);
}

double Score::run(Obj<RNA> rna) {
	_distAnal->readRNA(rna);
	_dihAnal->readRNA(rna);
	_distScore = _distAnal->scoring();
	_dihScore = _dihAnal->scoring();
	//cerr << rna->name << ' ' << _distScore << ' ' << _dihScore << endl;
	//double *score = _distAnal->getScore();
	//cerr << rna->name << ' ' << score[0] << ' ' << score[1] << ' ' << score[2] << ' ' << score[3] << ' ' << dihScore << endl;
	//cerr << _dist_bin << ' ' << _dih_bin << endl;
	//cerr << _distWeight << ' ' << _dihWeight << endl;
	//cerr << _distScore << ' ' << _dihScore << endl;
	_score = _constant + _distWeight * _distScore + _dihWeight * _dihScore;
	return _score;
}

void Score::trainDistPar(string filename) {
	string str;
	DistAnal distAnal(cutoff, _dist_bin);
	ifstream ifile(filename.c_str());
	int n = 0;
	while (ifile >> str) {
		n++;
		cerr << n << ". train: " << str << endl;
		Obj<RNA> rna(new RNA(str));
		distAnal.readRNA(rna);
		distAnal.train();
	}
	ifile.close();
	distAnal.printObsParm();
}

void Score::trainDihPar(string filename) {
	string str;
	DihAnal dihAnal(_dih_bin);
	ifstream ifile(filename.c_str());
	int n = 0;
	while (ifile >> str) {
		n++;
		cerr << n << ". train: " << str << endl;
		Obj<RNA> rna(new RNA(str));
		dihAnal.readRNA(rna);
		dihAnal.train();
	}
	ifile.close();
	dihAnal.printParm();
}

