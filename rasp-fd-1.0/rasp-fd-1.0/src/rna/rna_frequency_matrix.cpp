#include "rna_potential.h"

void usage(char *);

int main (int argc, char * const argv[]) {

	bool
		atom_type_file = false,
		distance_class_file = false,
		out_file = false;

    char
		*prgn = argv[0],
		*atom_type_filename,
		*distance_class_filename,
		*out_filename,
		*filename,
		chn = 'A';
	
	int
		k = 0;
	
	if (argc == 1) {
		usage(prgn);
		exit(1);
	}


    //Option processor
	while (argc > 1 && argv[1][0] == '-') {
		if (argv[1][1] == 'a') {
			atom_type_file = true;
			if (argv[1][2] == '\0') {
				atom_type_filename = argv[2];
				++argv;
				--argc;
			} else {
				atom_type_filename = &argv[1][2];
			}
		} else if (argv[1][1] == 'd') {
			distance_class_file = true;
			if (argv[1][2] == '\0') {
				distance_class_filename = argv[2];
				++argv;
				--argc;
			} else {
				distance_class_filename = &argv[1][2];
			}
		} else if (argv[1][1] == 'k') {
			if (argv[1][2] == '\0') {
				k = tnstring::string2int(argv[2]);
				++argv;
				--argc;
			} else {
				k = tnstring::string2int(&argv[1][2]);
			}
		} else if (argv[1][1] == 'p') {
			if (argv[1][2] == '\0') {
				filename = argv[2];
				++argv;
				--argc;
			} else {
				filename = &argv[1][2];
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
		} else {
			std::cerr << "Bad option \"-" << argv[1][1] << "\"" << std::endl;
			usage(prgn);
			exit(1);
		}
		++argv;
		--argc;
	}

	//std::cerr << argv[1] << std::endl;

	std::ifstream
		atom_type_ifstream,
		distance_class_ifstream,
		pdb_file;

	int
		n_groups = 0,
		n_types = 0,
		n_dist_classes = 0,
		n_pt,
		size,
		ind,
		x,y,z,
		numscan;
		
	float
		d,
		max_d;

	tnpdb::pdb
		*input;

	std::vector<tnpot::point_type>
		nu_point_type;
	
	tnpdb::point
		coords;

	matrix
		*mtx,
		*pot;
		    
	//char
	//	*filename = argv[1];

			
	std::vector<std::string>
		labels;
			
	std::ostringstream
		ostr;
	
	std::ofstream
		out;

	
	if (atom_type_file) {
		atom_type_ifstream.open(atom_type_filename, std::ifstream::in);
		if (atom_type_ifstream.fail() || !atom_type_ifstream.is_open()) {
			std::cerr << "ERROR: Cannot open file \"" << atom_type_filename << "\"." << std::endl;
			exit(1);
		}
		tnpot::load_atom_type_file(atom_type_ifstream);
		n_types = tnpot::number_of_atom_types_per_label(tnpdb::residue::NU);
		//std::cerr << n_types << std::endl;
	} else {
		std::cerr << "ERROR: You have not attached any atom type filename." << std::endl;
		exit(1); 
	}

	if (distance_class_file) {
		distance_class_ifstream.open(distance_class_filename, std::ifstream::in);
		if (distance_class_ifstream.fail() || !distance_class_ifstream.is_open()) {
			std::cerr << "ERROR: Cannot open file \"" << distance_class_filename << "\"." << std::endl;
			exit(1);
		}
		tnpot::load_distance_class_file(distance_class_ifstream);
		n_dist_classes = tnpot::number_of_distance_classes();
		max_d = tnpot::max_distance();
		//std::cerr << n_dist_classes << " " << max_d << std::endl;
	} else {
		std::cerr << "ERROR: You have not attached any distance class filename." << std::endl;
		exit(1);
	}


	ostr.str("");
	ostr << "+---------------------------------------------------";
	for (unsigned int i = 0; filename[i] != '\0'; ++i) ostr << "-";
	ostr << "---------------------------------+" << std::endl;
	std::cerr << ostr.str();
	std::cerr << "|                                 Processing File: " << filename << "                                  |" << std::endl;
	std::cerr << ostr.str();
	
	//Loading PDB file to memory
	pdb_file.open(filename, std::ifstream::in);
	if (pdb_file.fail() || !pdb_file.is_open()) {
		std::cerr << "ERROR: Cannot open file \"" << filename << "\"." << std::endl;
		exit(1);
	}
	input = new tnpdb::pdb(pdb_file);
	pdb_file.close();
		
	//Getting only valid (atom) types
	nu_point_type = tnpot::get_point_types(*input, tnpdb::residue::NU);
	//The number of (atom) types
	n_pt = nu_point_type.size();
	
	mtx = new matrix(4, k + 1, n_types, n_types, n_dist_classes);
	
	tnpot::get_rna_contact_matrix(mtx, &nu_point_type);

	if (out_file) {
		out.open(out_filename, std::iostream::out);
		tnpot::print_rna_contact_matrix(*mtx,out);
		out.close();
	} else {
		tnpot::print_rna_contact_matrix(*mtx, std::cout);
	}

    return 0;
}


void usage(char *prgn){
    std::cout <<
		"Usage: " << prgn << " -a <filename> -d <filename> -k <int> [-o <filename>] -p <filename>\n"
		"\tParameters:\n"
		"\t-a <filename>       File containing atom type definition.\n"
		"\t-d <filename>       File containing distance class definition.\n"
		"\t-k <int>            Interactions between residue i and residues >= i + k + 1.\n"
		"\t-p <filename>       PDB file containing RNA coordinates\n\n"
		"\tOptions:\n"
		"\t-o <filename>       File for writing output\n\n"
		"\tRNA FREQUENCIES\n"
		"\tThis programme calculates a frequency matrix of pairwise contacts in a PDB file\n"
		"\tcontaining RNA coordinates according to an atom type and a distance class definition.\n\n"
		"\tTomas Norambuena A. <tanoramb@puc.cl>\n";
	
	return;
}


