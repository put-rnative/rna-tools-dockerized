#ifndef PAR_H
#define PAR_H

#include "std.h"

class Par {
public:
	Par() {}
	Par(int argc, char **argv) {
		read(argc, argv);
	}
	std::vector<std::string> operator [](std::string key) {
		return _par[key];
	}
	int count(std::string key) {
		return _par.count(key);
	}
	void read(int, char **);
	map<std::string, std::vector<std::string>> _par;
};

#endif
