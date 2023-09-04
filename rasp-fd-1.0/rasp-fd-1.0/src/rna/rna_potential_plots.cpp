#include "rna_potential.h"

void usage(char *);

int main (int argc, char * const argv[]) {

	bool
		out_file = false;

    char
		*prgn = argv[0],
		*out_filename;
	
	if (argc == 1) {
		usage(prgn);
		exit(1);
	}


    //Option processor
	while (argc > 1 && argv[1][0] == '-') {
		if (argv[1][1] == 'o') {
			out_file = true;
			if (argv[1][2] == '\0') {
				out_filename = argv[2];
				++argv;
				--argc;
			} else {
				out_filename = &argv[1][2];
			}
		} else {
			std::cerr << "Bad option" << std::endl;
			usage(prgn);
			exit(1);
		}
		++argv;
		--argc;
	}

	//std::cerr << argv[1] << std::endl;

	std::ifstream
		potential_list_ifstream,
		potential_ifstream;

	int
		n = 0,
		numscan,
		k;

	std::vector<matrix>
		mtxs;
		    
	char
		*potential_list_filename = argv[1],
		potential_filename[200];
			
	std::ostringstream
		ostr;
	
	std::ofstream
		out;
		
	std::string
		line;
		
	std::vector<int>
		dim,
		ks;
		

	potential_list_ifstream.open(potential_list_filename, std::ifstream::in);
	
	while (getline(potential_list_ifstream, line, '\n')) {
		line = tnstring::trim(line);
		if (line[0] != '#' && line != "") {
				
			numscan = sscanf(line.c_str(), "%s\t%d", potential_filename, &k);
			
			ostr.str("");
			ostr << "+---------------------------------------------------";
			for (unsigned int i = 0; potential_filename[i] != '\0'; ++i) ostr << "-";
			ostr << "---------------------------------+" << std::endl;
			std::cerr << ostr.str();
			std::cerr << "|                                 Processing File: " << potential_filename << "                                  |" << std::endl;
			std::cerr << ostr.str();
				
			//Loading contact matrix file to memory
			potential_ifstream.open(potential_filename, std::ifstream::in);
			mtxs.push_back(tnpot::upload_rna_potential_matrix(potential_ifstream));
			potential_ifstream.close();
				
			ks.push_back(k);
		}
	}
	
	potential_list_ifstream.close();
	//for (unsigned int i = 0; i < ks.size(); i++) std::cerr << ks[i] << std::endl;	
	tnpot::asgl_potential_plots_by_k(&mtxs,ks);

    return 0;
}


void usage(char *prgn) {
    std::cout <<
		"Usage: " << prgn << " <filename>\n"
		"\tParameters:\n"
		"\t<filename>           File containing a list of RNA potential matrix paths with k, in the following\n"
		"\t                     format:\n\n"
		"\t                      [path_to_potential_matrix].  [k]\n\n"
		"\tPOTENTIAL PLOT BY K\n"
		"\tThis programme generates data for plotting potentials according to k. Plots are in TOP format (ASGL).\n\n"
		"\tTomas Norambuena A. <tanoramb@puc.cl>\n";
	
	return;
}


