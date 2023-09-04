#include "rna_potential.h"

void usage(char *);

int main (int argc, char * const argv[]) {

	bool
		atom_type_file = false,
		distance_class_file = false,
		potential_matrix_file = false,
		pdb_file = false,
		out_file = false;

    char
		*prgn = argv[0],
		*atom_type_filename,
		*distance_class_filename,
		*potential_matrix_filename,
		*out_filename,
		*pdb_filename,
		chn = 'A';
	
	int
		k = 0,
		z = 200;
	
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
		} else if (argv[1][1] == 'e') {
			potential_matrix_file = true;
			if (argv[1][2] == '\0') {
				potential_matrix_filename = argv[2];
				++argv;
				--argc;
			} else {
				potential_matrix_filename = &argv[1][2];
			}
		} else if (argv[1][1] == 'z') {
			if (argv[1][2] == '\0') {
				z = tnstring::string2int(argv[2]);
				++argv;
				--argc;
			} else {
				z = tnstring::string2int(&argv[1][2]);
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
	
	
	double
		nrg,
		nrg_norm,
		nrg_mean,
		nrg_sd,
		z_score;
		
	int
		n_contacts;
		
	tnpot::get_nrg_n_contacts(&nrg, &n_contacts, *input, &pot);
	tnpot::get_nrg_mean_nrg_sd(&nrg_mean,&nrg_sd,*input,&pot,z);

	nrg_norm = nrg / n_contacts;
	if (nrg_sd == 0.0) z_score = 0.0;
	else z_score = (nrg - nrg_mean) / nrg_sd;

	if (out_file) {
		out.open(out_filename, std::iostream::out);
		out << nrg << "\t" << n_contacts << "\t" << nrg_norm << "\t" << nrg_mean << "\t" << nrg_sd << "\t" << z_score <<std::endl;
		out.close();
	} else {
		std::cout << nrg << "\t" << n_contacts << "\t" << nrg_norm << "\t" << nrg_mean << "\t" << nrg_sd << "\t" << z_score << std::endl;
	}

    return 0;
}


void usage(char *prgn){
    std::cout <<
		"Usage: " << prgn << " -a <filename> -d <filename> -e <filename> [-z <int>] [-o <filename>] -p <filename>\n"
		"\tParameters:\n"
		"\t-a <filename>       File containing atom type definition.\n"
		"\t-d <filename>       File containing distance class definition.\n"
		"\t-e <filename>       File containing the corresponding potential matrix.\n"
		"\t-p <filename>       PDB file containing RNA coordinates\n\n"
		"\tOptions:\n"
		"\t-z <int>            Number of randomizations for calculating Z-Score. Default = 200. If 0 then Z-Score = 0.\n"
		"\t-o <filename>       File for writing output\n\n"
		"\tRNA POTENTIAL EVALUATOR\n"
		"\tThis programme calculates the energy for a PDB file, according to the potential matrix supplied:\n\n"
		"\t [PDB Energy]  [Number of Contacts]  [Normalized Energy]  [Mean Energy]  [SD Energy]  [Z-Score]\n\n"
		"\tTomas Norambuena A. <tanoramb@puc.cl>\n";
	
	return;
}


