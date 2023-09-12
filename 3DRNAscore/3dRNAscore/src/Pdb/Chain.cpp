#include "Chain.h"

Chain::Chain() {
	name = "A";
}

Chain::Chain(vector<string> &lines, string rnaName) {
	/* set name */
	name += lines[0][21];

	/* set rnaName */
	this->rnaName = rnaName;

	/* set residues */
	vector<string> strings;
	for (int i = 0; i < (int) lines.size(); i++) {
		if (strings.size() != 0 && lines[i].substr(22, 5) != strings.back().substr(22, 5)) {
			Residue *residue = new Residue(strings, rnaName);
			strings.clear();
			residues.push_back(*residue);
			delete residue;
			/*
			if (temp == "/A" || temp == "/U" || temp == "/G" || temp == "/C") {
				residues.push_back(*residue);
				delete residue;
			} else if ((residue->name == "A" && residue->getAmount() != 22 && residue->getAmount() != 19) ||
				(residue->name == "U" && residue->getAmount() != 20 && residue->getAmount() != 17) ||
				(residue->name == "C" && residue->getAmount() != 20 && residue->getAmount() != 17) ||
				(residue->name == "G" && residue->getAmount() != 23 && residue->getAmount() != 20)) {
				delete residue;
			} else {
				residues.push_back(*residue);
				delete residue;
			}
			*/
		}
		strings.push_back(lines[i]);
	}
	Residue *residue = new Residue(strings, rnaName);
	residues.push_back(*residue);
	delete residue;
	/*
	if (temp == "/A" || temp == "/U" || temp == "/G" || temp == "/C") {
		residues.push_back(*residue);
		delete residue;
	} else if ((residue->name == "A" && residue->getAmount() != 22 && residue->getAmount() != 19) ||
		(residue->name == "U" && residue->getAmount() != 20 && residue->getAmount() != 17) ||
		(residue->name == "C" && residue->getAmount() != 20 && residue->getAmount() != 17) ||
		(residue->name == "G" && residue->getAmount() != 23 && residue->getAmount() != 20)) {
		delete residue;
	} else {
		residues.push_back(*residue);
		delete residue;
	}
	*/
}

void Chain::push(Residue *residue) {
	residues.push_back(*residue);
}

void Chain::push(Residue &residue) {
	residues.push_back(residue);
}

Residue &Chain::operator [](int n) {
	return residues[n];
}

