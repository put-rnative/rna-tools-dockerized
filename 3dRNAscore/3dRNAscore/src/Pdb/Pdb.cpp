#include "Pdb.h"

Pdb::Pdb() {
	name = "";
}

Pdb::Pdb(char *pdbfile) {
	string str(pdbfile);
	readPDB(pdbfile);
}

Pdb::Pdb(string pdbfile) {
	readPDB(pdbfile);
}

Pdb *Pdb::copy() {
	Pdb *pdb = new Pdb;
	pdb->name = name;
	pdb->len = len;
	pdb->chains = chains;
	return pdb;
}

void Pdb::readPDB(string pdbfile) {
	/* set name */
	if (pdbfile.size() <= 4 || pdbfile.substr(pdbfile.size() - 4, 4).compare(".pdb")) {
		cerr << "Please give me a file ended with '.pdb'\n";
		exit(1);
	}
	name = pdbfile.substr(0, pdbfile.size() - 4);
	/*
	string tmp;
	for (int i = name.size() - 1; i >= 0 && name[i] != '\\' && name[i] != '/'; i--) {
		tmp += name[i];
	}
	name = "";
	for (int i = tmp.size() - 1; i >= 0; i--) {
		name += tmp[i];
	}*/

	/* set chains */
	ifstream ifile(pdbfile.c_str());
	if (!ifile) {
		cout << "open file \"" << pdbfile << "\" failed!" << endl;
		exit(1);
	}
	string line;
	vector<string> lines;
	int n = 0;
	while (getline(ifile, line, '\n')) {
		if (!line.compare(0, 4, "ATOM")) {
			n++;
			if (lines.size() != 0 && line[21] != lines.back()[21]) {
				Obj<Chain> chain = new Chain(lines, name);
				lines.clear();
				chains.push_back(*chain);
			}
			lines.push_back(line);
		}
	}

	if (n == 0) {
		cerr << "The file '" << pdbfile << "' has nothing!" << endl;
		exit(1);
	}

	Obj<Chain> chain = new Chain(lines, name);
	chains.push_back(*chain);
	ifile.close();

	/* set length */
	setLen();
}

Chain &Pdb::operator [](int n) {
	return chains[n];
}

void Pdb::push(Chain *chain) {
	chains.push_back(*chain);
}

void Pdb::push(Chain &chain) {
	chains.push_back(chain);
}

string Pdb::getSeq() {
	string str;
	for (int i = 0; i < (int) chains.size(); i++) {
		for (int j = 0; j < (int) chains[i].residues.size(); j++) {
			string temp = chains[i].residues[j].name;
			for (int k = 0; k < (int) temp.size(); k++) {
				if (temp[k] == 'A' || temp[k] == 'U' || temp[k] == 'C' || temp[k] == 'G') {
					str += temp[k];
					break;
				}
			}
		}
	}
	return str;
}

void Pdb::print() {
	int i, j, k, flag, flag2;
	double x, y, z;

	cout << fixed << setprecision(3);
	for (i = 0, flag = 0, flag2 = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			flag2++;
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				if (j == 0) {
					name = chains[i].residues[j].atoms[k].name;
					if (name == "P" || name == "O1P" || name == "O2P") {
						continue;
					}
				}
				flag++;
				x = chains[i].residues[j].atoms[k].x;
				y = chains[i].residues[j].atoms[k].y;
				z = chains[i].residues[j].atoms[k].z;
				cout << "ATOM  " << setw(5) << flag << "  " << setw(4) << left << chains[i].residues[j].atoms[k].name << right << setw(3) << chains[i].residues[j].name << setw(2) << char('A' + i) << setw(4) << flag2 << setw(12) << x << setw(8) << y << setw(8) << z << '\n';
			}
		}
		cout << "TER\n";
	}
}

void Pdb::printAsDNA() {
	int i, j, k, flag, flag2;
	double x, y, z;

	cout << fixed << setprecision(3);
	for (i = 0, flag = 0, flag2 = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			flag2++;
			if (chains[i].residues[j].name == "A") {
				chains[i].residues[j].name = "DA";
			} else if (chains[i].residues[j].name == "U") {
				chains[i].residues[j].name = "DT";
			} else if (chains[i].residues[j].name == "G") {
				chains[i].residues[j].name = "DG";
			} else if (chains[i].residues[j].name == "C") {
				chains[i].residues[j].name = "DC";
			}
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				if (j == 0) {
					name = chains[i].residues[j].atoms[k].name;
					if (name == "P" || name == "O1P" || name == "O2P") {
						continue;
					}
				}
				if (chains[i].residues[j].atoms[k].name == "O2*") continue;
				flag++;
				x = chains[i].residues[j].atoms[k].x;
				y = chains[i].residues[j].atoms[k].y;
				z = chains[i].residues[j].atoms[k].z;
				cout << "ATOM  " << setw(5) << flag << "  " << setw(4) << left << chains[i].residues[j].atoms[k].name << right << setw(3) << chains[i].residues[j].name << setw(2) << char('A' + i) << setw(4) << flag2 << setw(12) << x << setw(8) << y << setw(8) << z << '\n';
			}
		}
		cout << "TER\n";
	}
}

void Pdb::write(string pdbname) {
	int i, j, k, flag, flag2;
	ofstream ofile;
	double x, y, z;

	ofile.open(pdbname.c_str());
	ofile << fixed << setprecision(3);
	for (i = 0, flag = 0, flag2 = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			flag2++;
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				if (j == 0) {
					name = chains[i].residues[j].atoms[k].name;
					if (name == "P" || name == "O1P" || name == "O2P") {
						continue;
					}
				}
				flag++;
				x = chains[i].residues[j].atoms[k].x;
				y = chains[i].residues[j].atoms[k].y;
				z = chains[i].residues[j].atoms[k].z;
				ofile << "ATOM  " << setw(5) << flag << "  " << setw(4) << left << chains[i].residues[j].atoms[k].name << right << setw(3) << chains[i].residues[j].name << setw(2) << char('A' + i) << setw(4) << flag2 << setw(12) << x << setw(8) << y << setw(8) << z << '\n';
			}
		}
		ofile << "TER\n";
	}
	ofile.close();
}

void Pdb::setLen() {
	int i, j, flag;

	for (flag = 0, i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			flag++;
		}
	}
	len = flag;
}

inline int Pdb::getLen() {
//	return len;
	
	int i, j, flag;

	for (flag = 0, i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			flag++;
		}
	}

	return flag;
	
}

void Pdb::updateChains(string ss) {
	// compute counts of chain
	int temp = 0;
	for (int i = 0; i < (int) ss.size(); i++) {
		if (ss[i] == '&' || (ss[i] == ')' && i - 1 >= 0 && ss[i - 1] == '(')) {
			temp++;
		}
	}
	temp++;

	// construct new chains
	Chain *chain = new Chain[temp];
	int i1 = 0;
	for (int j = 0, k = 0; i1 < (int) chains.size(); i1++) {
		for (int i2 = 0; i2 < (int) chains[i1].residues.size();) {
			if (ss[k] == '&' || (ss[k] == ')' && k - 1 >= 0 && ss[k - 1] == '(')) {
				k++;
				j++;
				continue;
			}
			chain[j].residues.push_back(chains[i1].residues[i2]);
			k++;
			i2++;
		}
	}
	chains.clear();
	for (int i1 = 0; i1 < temp; i1++) {
		chains.push_back(chain[i1]);
	}
	delete [] chain;
}

int Pdb::totalAtoms() {
	int i, j, k, flag;

	for (flag = 0, i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
				flag++;
			}
		}
	}

	return flag;
}

double Pdb::getDist(int a, int b) {
	int i, j, k, flag;
	double x1, y1, z1, x2, y2, z2;

	if (a > b || b > getLen() || a > getLen()) {
		return -1;
	}

	for (flag = 0, i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			flag++;
			if (flag == a) {
				for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
					if (chains[i].residues[j].atoms[k].name == "O5*") {
						x1 = chains[i].residues[j].atoms[k].x;
						y1 = chains[i].residues[j].atoms[k].y;
						z1 = chains[i].residues[j].atoms[k].z;
					}
				}
			}
			if (flag == b) {
				for (k = 0; k < (int) chains[i].residues[j].atoms.size(); k++) {
					if (chains[i].residues[j].atoms[k].name == "O3*") {
						x2 = chains[i].residues[j].atoms[k].x;
						y2 = chains[i].residues[j].atoms[k].y;
						z2 = chains[i].residues[j].atoms[k].z;
					}
				}
			}
		}
	}

	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

void Pdb::setResNum() {
	int i, j;
	int num;

	for (num = 0, i = 0; i < (int) chains.size(); i++) {
		for (j = 0; j < (int) chains[i].residues.size(); j++) {
			num++;
			chains[i].residues[j].num = num;
		}
	}
}

Pdbs::Pdbs() {
}

Pdbs::~Pdbs() {
	for (int i = 0; i < (int) PdbList.size(); i++) {
		delete PdbList[i];
	}
}

int Pdbs::getLen() {
	return PdbList.size();
}

void Pdbs::resize(int n) {
	PdbList.resize(n);
}

Pdb *Pdbs::at(int n) {
	return PdbList[n];
}

Pdb *Pdbs::operator [](int n) {
	return PdbList[n];
}

void Pdbs::push(Pdb *pdb) {
	PdbList.push_back(pdb);
}



