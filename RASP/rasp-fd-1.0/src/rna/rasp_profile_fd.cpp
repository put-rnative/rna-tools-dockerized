#include "rna_potential.h"

void usage(char *);

int main (int argc, char * const argv[]) {

	bool
		atom_type_file = true,
		distance_class_file = true,
		potential_matrix_file = true,
		user_potential = false,
		pdb_file = false,
		out_file = false,
		rasmol_scr = false;

	char
		*prgn = argv[0],
		atom_type_filename[1000],
		distance_class_filename[1000],
		potential_matrix_filename[1000],
		*out_filename,
		*pdb_filename,
		*rasp_c,
		chn = 'A';
	
	
	int
		k = 0,
		win = 1;
		
	float
		rasmol_cutoff = 0.0;
	
	if (argc == 1) {
		usage(prgn);
		exit(1);
	}


    //Option processor
	while (argc > 1 && argv[1][0] == '-') {
		if (argv[1][1] == 'e') {
			user_potential = true;
			if (argv[1][2] == '\0') {
				rasp_c = argv[2];
				++argv;
				--argc;
			} else {
				rasp_c = &argv[1][2];
			}
		} else if (argv[1][1] == 'w') {
			if (argv[1][2] == '\0') {
				win = tnstring::string2int(argv[2]);
				++argv;
				--argc;
			} else {
				win = tnstring::string2int(&argv[1][2]);
			}
		} else if (argv[1][1] == 'o') {
			out_file = true;
			if (argv[1][2] == '\0') {
				out_filename = argv[2];
				++argv;
				--argc;
			} else {
				out_filename = &argv[1][2];
			}
		} else if (argv[1][1] == 'p') {
			pdb_file = true;
			if (argv[1][2] == '\0') {
				pdb_filename = argv[2];
				++argv;
				--argc;
			} else {
				pdb_filename = &argv[1][2];
			}
		} else if (argv[1][1] == 'r') {
			rasmol_scr = true;
			/*if (argv[1][2] == '\0') {
				rasmol_cutoff = tnstring::string2float(argv[2]);
				++argv;
				--argc;
			} else {
				rasmol_cutoff = tnstring::string2float(&argv[1][2]);
			}
			*/
		} else {
			std::cerr << "Bad option \"-" << argv[1][1] << "\"" << std::endl;
			usage(prgn);
			exit(1);
		}
		++argv;
		--argc;
	}

	//std::cerr << argv[1] << std::endl;
	if (!user_potential) {
		rasp_c = new char[4];
		strcpy(rasp_c,"all");
	}

	if (getenv("RASP") == NULL) {
		std::cerr << "You must set the environment variable 'RASP' according to the instructions in the README.txt file." << std::endl;
		exit(1);
	}
	
	char
		*env_c = getenv("RASP");

	std::string
		rasp(rasp_c),
		env(env_c),
		path;

	path = env + "/lib/distances.dat";
	strcpy(distance_class_filename,path.c_str());
	
	if (tnstring::toLower(rasp) == "all") {
		path = env + "/lib/types_all.dat";
		strcpy(atom_type_filename, path.c_str());
		path = env + "/lib/all.nrg";
		strcpy(potential_matrix_filename, path.c_str());
	} else if (tnstring::toLower(rasp) == "bbr") {
		path = env + "/lib/types_bbr.dat";
		strcpy(atom_type_filename,path.c_str());
		path = env + "/lib/bbr.nrg";
		strcpy(potential_matrix_filename,path.c_str());
	} else if (tnstring::toLower(rasp) == "bb") {
		path = env + "/lib/types_bb.dat";
		strcpy(atom_type_filename,path.c_str());
		path = env + "/lib/bb.nrg";
		strcpy(potential_matrix_filename,path.c_str());
	} else if (tnstring::toLower(rasp) == "c3") {
		path = env + "/lib/types_c3.dat";
		strcpy(atom_type_filename,path.c_str());
		path = env + "/lib/c3.nrg";
		strcpy(potential_matrix_filename,path.c_str());
	} else {
		std::cerr << "Bad potential" << std::endl;
		usage(prgn);
		exit(1);
	}
	
	if (win % 2 == 0 || win > 9 || win < 1) {
		std::cerr << "ERROR: Bad window size" << std::endl;
		usage(prgn);
		exit(1);
	}

	std::ifstream
		atom_type_ifstream,
		distance_class_ifstream,
		potential_matrix_ifstream,
		pdb_ifstream;

	tnpdb::pdb
		*input;

	matrix
		pot;
			
	std::vector<std::string>
		labels;
			
	std::ostringstream
		ostr;
	
	std::ofstream
		out;

	
	if (atom_type_file) {
		atom_type_ifstream.open(atom_type_filename, std::ifstream::in);
		if (!atom_type_ifstream.is_open() || atom_type_ifstream.bad()) {
			std::cerr << "ERROR: Cannot open file " << atom_type_filename << std::endl;
			exit(1);
		}
		tnpot::load_atom_type_file(atom_type_ifstream);
		//n_types = tnpot::number_of_atom_types_per_label(tnpdb::residue::NU);
		atom_type_ifstream.close();
	} else {
		std::cerr << "ERROR: You have not attached any atom type filename." << std::endl;
		exit(1); 
	}

	if (distance_class_file) {
		distance_class_ifstream.open(distance_class_filename, std::ifstream::in);
		if (!distance_class_ifstream.is_open() || distance_class_ifstream.bad()) {
			std::cerr << "ERROR: Cannot open file " << distance_class_filename << std::endl;
			exit(1);
		}
		tnpot::load_distance_class_file(distance_class_ifstream);
		//n_dist_classes = tnpot::number_of_distance_classes();
		//max_d = tnpot::max_distance();
		distance_class_ifstream.close();
	} else {
		std::cerr << "ERROR: You have not attached any distance class filename." << std::endl;
		exit(1);
	}

	if (potential_matrix_file) {
		potential_matrix_ifstream.open(potential_matrix_filename, std::ifstream::in);
		if (!potential_matrix_ifstream.is_open() || potential_matrix_ifstream.bad()) {
			std::cerr << "ERROR: Cannot open file " << potential_matrix_filename << std::endl;
			exit(1);
		}
		pot = tnpot::upload_rna_potential_matrix(potential_matrix_ifstream);
		potential_matrix_ifstream.close();
	} else {
		std::cerr << "ERROR: You have not attached any potential matrix filename." << std::endl;
		exit(1);
	}
	
	if (pdb_file) {

		ostr.str("");
		ostr << "+---------------------------------------------------";
		for (unsigned int i = 0; pdb_filename[i] != '\0'; ++i) ostr << "-";
		ostr << "---------------------------------+" << std::endl;
		std::cerr << ostr.str();
		std::cerr << "|                                 Processing File: " << pdb_filename << "                                  |" << std::endl;
		std::cerr << ostr.str();
	
		//Loading PDB file to memory
		pdb_ifstream.open(pdb_filename, std::ifstream::in);
			if (!pdb_ifstream.is_open() || pdb_ifstream.bad()) {
			std::cerr << "ERROR: Cannot open file " << pdb_filename << std::endl;
			exit(1);
		}
		input = new tnpdb::pdb(pdb_ifstream);
		pdb_ifstream.close();
	} else {
		std::cerr << "ERROR: You have not supplied any PDB filename." << std::endl;
		exit(1);
	}
	
		
	int
		size,
		i;
		
	std::vector<tnpot::pick>
		profile = tnpot::get_nrg_profile(*input, &pot, win);

	size = profile.size();

	if (out_file) {
		out.open(out_filename, std::iostream::out);
		for (i = 0; i < size; i++) out << tnstring::trim(profile[i]._resname) << "\t" << profile[i]._resindex << "\t" << profile[i]._chn << "\t" << profile[i]._value << std::endl;
		out.close();
	} else {
		for (i = 0; i < size; i++) std::cout << tnstring::trim(profile[i]._resname) << "\t" << profile[i]._resindex << "\t" << profile[i]._chn << "\t" << profile[i]._value << std::endl;
	}
	
	if (rasmol_scr) {
		
		print_rasmol_script(&profile);
		
	}

    return 0;
}


void usage(char *prgn){
	char ESC=27;
    std::cout <<
		"Usage: " << prgn << " -e <string> [-w <int>] [-o <filename>] [-r] -p <filename>\n"
		"\tParameters:\n"
		"\t-e <string>         The potential to use: c3, bb, bbr, all. Default: all.\n"
		"\t-p <filename>       PDB file containing RNA coordinates.\n\n"
		"\tOptions:\n"
		"\t-w <int>            Window to average profile score. Possible values: 1, 3, 5, 7, 9. Default: 1\n"
		"\t-r                  Print a script for Rasmol visualization (called 'profile.scr'). You must execute it this way:\n\n"
		"\t                      rasmol -script profile.scr path_to_pdb_file\n\n"
		"\t-o <filename>       File for writing output\n\n"
		"\t" << ESC <<"[1mRNA POTENTIAL PROFILE" << ESC << "[0m\n"
		"\tThis programme calculates the energy profile for a RNA PDB file, according to the potential selected. The output\n"
		"\tis as follow:\n\n"
		"\t [Residue Name]  [Index]  [Chain ID]  [Energy]\n\n"
		"\t" << ESC <<"[1mTomas Norambuena A." << ESC << "[0m  <tanoramb@puc.cl>\n"
		"\tMolecular Bioinformatics Laboratory\n"
		"\tPontificia Universidad Catolica de Chile\n";
		
	
	return;
}


