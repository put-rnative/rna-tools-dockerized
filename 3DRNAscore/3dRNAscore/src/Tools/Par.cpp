#include "Par.h"

void Par::read(int argc, char **argv) {
	/*
	for (int i = 0; i < argc; i++) {
		std::smatch sm;
		if (std::regex_match(string(argv[i]), sm, std::regex("^-(\\w+):(\\w+)$"))) {
			_par[sm[1]] = sm[2];
		} else if (std::regex_match(string(argv[i]), sm, std::regex("^-(\\w+)$"))) {
			_par[sm[1]] = "none";
		}
	}
	*/
	std::string key;
	std::vector<std::string> values;
	int n = 0;
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (n != 0) {
				_par[key] = values;
			}
			values.clear();
			n++;

			std::string str(argv[i]);
			if (strlen(argv[i]) == 1) {
				std::cerr << "Par::read error! The parameter should not be '-'. " << std::endl;
				exit(1);
			} else {
				key = str.substr(1, str.size() - 1);
			}
		} else {
			values.push_back(argv[i]);
		}
	}
	if (n != 0) {
		_par[key] = values;
	}
}




