#include "../Score.h"

int main(int argc, char **argv) {
	Par par(argc, argv);
	if (par.count("-score:list")) {
		string str;
		Obj<Score> score;
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(par["-score:list"][1]);
		}

		ifstream ifile(par["-score:list"][0]);
		while (ifile >> str) {
			Obj<RNA> rna = new RNA(str);
			cout << str << ' ' << score->run(rna) << endl;
		}
		ifile.close();
	}
	if (par.count("s:l")) {
		string str;
		Obj<Score> score;
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(par["s:l"][1]);
		}

		ifstream ifile(par["s:l"][0]);
		while (ifile >> str) {
			Obj<RNA> rna = new RNA(str);
			cout << str << ' ' << score->run(rna) << endl;
		}
		ifile.close();
	}
	if (par.count("-score")) {
		Obj<Score> score;
		Obj<RNA> rna = new RNA(par["-score"][0]);
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(par["-score"][1]);
		}
		cout << score->run(rna) << endl;
	}
	if (par.count("s")) {
		Obj<Score> score;
		Obj<RNA> rna = new RNA(par["s"][0]);
		if (argc == 3) {
			score = new Score();
		} else if (argc == 4) {
			score = new Score(par["s"][1]);
		}
		cout << score->run(rna) << endl;
	}
	if (par.count("h") || par.count("-help")) {
		std::cout << "Version: 3dRNAscore-1.0" << std::endl;
		std::cout << std::endl;
		std::cout << "Caution: Please first set the 'RNAscore' environment to the path 3dRNAscore directory." << std::endl;
		std::cout << std::endl;
		std::cout << "         It's better to check if the pdb files lack some atoms. We have provided a perl " << std::endl;
		std::cout << "         script named 'format.pl' in the 'lib' directory to help to check missing atoms. " << std::endl;
		std::cout << "         Just type './format.pl -check filename'. " << std::endl;
		std::cout << std::endl;
		std::cout << "Usage: 1   ./3dRNAscore {-s}|{--score} <pdb_file> [<parameter_file>]" << std::endl;
		std::cout << "               # score a single RNA named '<pdb_file>'." << std::endl;
		std::cout << "                 We've provided a sample parameter file in the 'lib' directory. The user " << std::endl;
		std::cout << "                 can control the parameters through the parameter file. " << std::endl;
		std::cout << "       2   ./3dRNAscore {-s:l}|{--score:list} <list_file> [<parameter_file>]" << std::endl;
		std::cout << "               # score a list of RNA whose name were listed in the file '<parameter_file>'" << std::endl;
		std::cout << "       3   ./3dRNAscore {-h}|{--help}" << std::endl;
		std::cout << "               # help information" << std::endl;
		std::cout << "       4   ./3dRNAscore {-v}|{--version}" << std::endl;
		std::cout << "               # version information" << std::endl;
		std::cout << std::endl;
		std::cout << "Please feel free to contact the 3dRNAscore team if you have any question!" << std::endl;
		std::cout << "wj_hust08@hust.edu.cn yjzhao.wh@gmail.com" << std::endl;
		std::cout << std::endl;
	}
	if (par.count("v") || par.count("-version")) {
		std::cout << "Version: 3dRNAscore-1.0" << std::endl;
		std::cout << std::endl;
	}
	return 0;
}
