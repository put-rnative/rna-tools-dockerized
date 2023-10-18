#include "DistAnal.h"

DistAnal::DistAnal(int cutoff, double interval, string reference_state) {
	num = NULL;
	type = NULL;
	ntLen = NULL;
	list = NULL;
	len = 0;

	obsProb = NULL;
	refProb = NULL;
	score = 0;

	this->interval = interval;
	this->cutoff = cutoff;
	this->reference_state = reference_state;
	bins = int(ceil(cutoff / interval));
	penalty = 0;

	obsParm = new int[85 * 85 * bins];
	obsProb = new double[85 * 85 * bins];
	for (int i = 0; i < 85 * 85 * bins; i++) {
		obsParm[i] = 0;
		obsProb[i] = 0;
	}
	refParm = new int[85 * 85 * bins];
	refProb = new double[85 * 85 * bins];
	for (int i = 0; i < 85 * 85 * bins; i++) {
		refParm[i] = 0;
		refProb[i] = 0;
	}
}

DistAnal::~DistAnal() {
	delete [] num;
	delete [] type;
	delete [] ntLen;
	for (int i = 0; i < len; i++) {
		delete [] list[i];
	}
	delete [] list;

	delete [] obsParm;
	delete [] refParm;

	delete [] obsProb;
	delete [] refProb;
}

void DistAnal::readRNA(Obj<RNA> rna) {
	delete [] num;
	delete [] type;
	delete [] ntLen;
	for (int i = 0; i < len; i++) {
		delete [] list[i];
	}
	delete [] list;

	len = rna->getLen();
	num = new int[len];
	type = new int[len];
	ntLen = new int[len];
	list = new Point *[len];

	Point *p = new Point;
	for (int i = 0, m = 0, n = 0; i < (int) rna->chains.size(); i++) {
		for (int j = 0; j < (int) rna->chains[i].residues.size(); j++, n++) {
			m++;
			ntLen[n] = (int) rna->chains[i].residues[j].atoms.size();
			string name = rna->chains[i].residues[j].name;
			if (name == "A") {
				type[n] = 1;
			} else if (name == "U") {
				type[n] = 2;
			} else if (name == "G") {
				type[n] = 3;
			} else if (name == "C") {
				type[n] = 4;
			}
			list[n] = new Point[ntLen[n]];
			for (int k = 0; k < (int) rna->chains[i].residues[j].atoms.size(); k++) {
				list[n][k].x = rna->chains[i].residues[j].atoms[k].x;
				list[n][k].y = rna->chains[i].residues[j].atoms[k].y;
				list[n][k].z = rna->chains[i].residues[j].atoms[k].z;
				if (rna->chains[i].residues[j].atoms[0].name != "P") {
					list[n][k].type = 3;
				} else {
					list[n][k].type = 0;
				}
				if (type[n] == 1) {
					list[n][k].type += k;
				} else if (type[n] == 2) {
					list[n][k].type += 22 + k;
				} else if (type[n] == 3) {
					list[n][k].type += 42 + k;
				} else {
					list[n][k].type += 65 + k;
				}
				if (rna->chains[i].residues[j].atoms[k].name == "O5*" && n != 0) {
					double dx = list[n][k].x - p->x;
					double dy = list[n][k].y - p->y;
					double dz = list[n][k].z - p->z;
					double dist = sqrt(dx * dx + dy * dy + dz * dz);
					if (dist > 4) {
						m += 10000;
					}
				} else if (rna->chains[i].residues[j].atoms[k].name == "O3*") {
					p->x = list[n][k].x;
					p->y = list[n][k].y;
					p->z = list[n][k].z;
				}
			}
			num[n] = m;
		}
	}
}

void DistAnal::train() {
	for (int i = 0; i < len; i++) {
		for (int j = i + 1; j < len; j++) {
			for (int k = 0; k < ntLen[i]; k++) {
				for (int l = 0; l < ntLen[j]; l++) {
					int type1 = list[i][k].type;
					int type2 = list[j][l].type;
					if (num[j] - num[i] == 1 && (((type1 >= 0 && type1 <= 11) || 
							(type1 >= 22 && type1 <= 33) || (type1 >= 42 && type1 <= 53) || 
							(type1 >= 62 && type1 <= 73)) && ((type2 >= 0 && type2 <= 11) || 
							(type2 >= 22 && type2 <= 33) || (type2 >= 42 && type2 <= 53) || 
							(type2 >= 62 && type2 <= 73)))) 
						continue;
					double temp = list[i][k].dist(&(list[j][l]));
					if (temp >= cutoff) continue;
					obsParm[(list[i][k].type * 85 + list[j][l].type) * bins + int(temp / interval)]++;
					obsParm[(list[j][l].type * 85 + list[i][k].type) * bins + int(temp / interval)]++;
					//refParm[int(temp / interval)] += 2;
				}
			}
		}
	}
}

void DistAnal::readObsParm(string filename) {
	ifstream ifile(filename.c_str());
	int temp;
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			ifile >> temp;
			if (!ifile) {
				break;
			}
			obsParm[i * bins + j] += temp;
		}
	}
	ifile.close();
	initObsProb();
}

void DistAnal::readRefParm(string filename) {
	if (reference_state == "average") {
		int *single_row_sum = new int[85 * 85];
		int *single_column_sum = new int[bins];
		int total = 0;
		for (int i = 0; i < 85 * 85; i++) {
			single_row_sum[i] = 0;
			for (int j = 0; j < bins; j++) {
				if (i == 0) {
					single_column_sum[j] = 0;
				}
				single_row_sum[i] += obsParm[i * bins + j];
				single_column_sum[j] += obsParm[i * bins + j];
				total += obsParm[i * bins + j];
			}
		}
		for (int i = 0; i < 85 * 85; i++) {
			for (int j = 0; j < bins; j++) {
				refParm[i * bins + j] = int(single_row_sum[i] * double(single_column_sum[j]) / total);
			}
		}
		delete [] single_row_sum;
		delete [] single_column_sum;
	} else {
		ifstream ifile(filename.c_str());
		int temp;
		for (int i = 0; i < 85 * 85; i++) {
			for (int j = 0; j < bins; j++) {
				ifile >> temp;
				refParm[i * bins + j] += temp;
			}
		}
		ifile.close();
	}
	initRefProb();
}

double DistAnal::scoring() {
	score = 0;
	for (int i = 0; i < len; i++) {
		for (int j = i + 1; j < len; j++) {
			for (int k = 0; k < ntLen[i]; k++) {
				for (int l = 0; l < ntLen[j]; l++) {
					int type1 = list[i][k].type;
					int type2 = list[j][l].type;
					if (num[j] - num[i] == 1 && (((type1 >= 0 && type1 <= 11) || (type1 >= 22 && type1 <= 33) || (type1 >= 42 && type1 <= 53) || (type1 >= 62 && type1 <= 73)) && ((type2 >= 0 && type2 <= 11) || (type2 >= 22 && type2 <= 33) || (type2 >= 42 && type2 <= 53) || (type2 >= 62 && type2 <= 73)))) continue;
					double temp = list[i][k].dist(&(list[j][l]));
					if (temp >= cutoff) continue;
					double a = obsProb[(list[i][k].type * 85 + list[j][l].type) * bins + int(temp / interval)];
					double b = obsProb[(list[j][l].type * 85 + list[i][k].type) * bins + int(temp / interval)];
					double c = refProb[(list[i][k].type * 85 + list[j][l].type) * bins + int(temp / interval)];
					double d = refProb[(list[j][l].type * 85 + list[i][k].type) * bins + int(temp / interval)];
					// double c = refProb[int(temp / interval)];
					if (a != 0 && c != 0) {
						if (((type1 > 11 && type1 < 22) || (type1 > 33 && type1 < 42) || (type1 > 54 && type1 < 62) || (type1 > 74 && type1 < 85)) && ((type2 > 11 && type2 < 22) || (type2 > 33 && type2 < 42) || (type2 > 54 && type2 < 62) || (type2 > 74 && type2 < 85))) {
							score -= 2.5 * log(a / c);
						} else {
							score -= log(a / c);
						}
					} else {
						score += penalty;
					}
					if (b != 0 && d != 0) {
						if (((type1 > 11 && type1 < 22) || (type1 > 33 && type1 < 42) || (type1 > 54 && type1 < 62) || (type1 > 74 && type1 < 85)) && ((type2 > 11 && type2 < 22) || (type2 > 33 && type2 < 42) || (type2 > 54 && type2 < 62) || (type2 > 74 && type2 < 85))) {
							score -= 2.5 * log(b / d);
						} else {
							score -= log(b / d);
						}
					} else {
						score += penalty;
					}
				}
			}
		}
	}
	score = score / (len * (len - 1));
	return score;
}

void DistAnal::initObsProb() {
	for (int i = 0; i < 85 * 85; i++) {
		double line_sum = 0;
		for (int j = 0; j < bins; j++) {
			line_sum += obsParm[i * bins + j];
		}
		for (int j = 0; j < bins; j++) {
			obsProb[i * bins + j] = double(obsParm[i * bins + j]) / line_sum;
		}
	}
}

void DistAnal::initRefProb() {
	for (int i = 0; i < 85 * 85; i++) {
		double line_sum = 0;
		for (int j = 0; j < bins; j++) {
			line_sum += refParm[i * bins + j];
		}
		for (int j = 0; j < bins; j++) {
			refProb[i * bins + j] = double(refParm[i * bins + j]) / line_sum;
		}
	}
}

void DistAnal::printObsParm() {
	int sum[bins];
	for (int i = 0; i < bins; i++) {
		sum[i] = 0;
	}
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			sum[j] += obsParm[i * bins + j];
			cout << obsParm[i * bins + j] << ' ';
		}
		cout << endl;
	}
	for (int i = 0; i < bins; i++) {
		cout << sum[i] << ' ';
	}
	cout << endl;
}

void DistAnal::printRefParm() {
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			cout << refParm[i * bins + j] << ' ';
		}
		cout << endl;
	}
}

void DistAnal::printObsProb() {
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			cout << obsProb[i * bins + j] << ' ';
		}
		cout << endl;
	}
}

void DistAnal::printRefProb() {
	for (int i = 0; i < 85 * 85; i++) {
		for (int j = 0; j < bins; j++) {
			cout << refProb[i * bins + j] << ' ';
		}
		cout << endl;
	}
}

double DistAnal::getScore() {
	return score;
}
