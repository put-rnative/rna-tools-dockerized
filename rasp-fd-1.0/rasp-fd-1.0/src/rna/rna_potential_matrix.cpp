#include "rna_potential.h"

void usage(char *);

int main (int argc, char * const argv[]) {

	bool
		print_asgl_plot = false,
		out_file = false,
		t_contact_matrix_file = false;

    char
		*prgn = argv[0],
		*contact_matrix_list_filename,
		*t_contact_matrix_filename,
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
		} else if (argv[1][1] == 'f') {
			t_contact_matrix_file = true;
			if (argv[1][2] == '\0') {
				t_contact_matrix_filename = argv[2];
				++argv;
				--argc;
			} else {
				t_contact_matrix_filename = &argv[1][2];
			}
		} else if (argv[1][1] == 'i') {
			if (argv[1][2] == '\0') {
				contact_matrix_list_filename = argv[2];
				++argv;
				--argc;
			} else {
				contact_matrix_list_filename = &argv[1][2];
			}
		} else if (argv[1][1] == 'g') {
			print_asgl_plot = true;
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
		contact_matrix_list_ifstream,
		contact_matrix_ifstream;

	int
		n = 0,
		numscan;

	matrix
		mtx,
		*sum,
		*pot;
		    
	char
		contact_matrix_filename[200];
			
	std::ostringstream
		ostr;
	
	std::ofstream
		out;
		
	std::string
		line;
		
	std::vector<int>
		dim;

	contact_matrix_list_ifstream.open(contact_matrix_list_filename, std::ifstream::in);
	if (contact_matrix_list_ifstream.fail() || !contact_matrix_list_ifstream.is_open()) {
		std::cerr << "ERROR: Cannot open file \"" << contact_matrix_list_filename << "\"." << std::endl;
		exit(1);
	}

	while (getline(contact_matrix_list_ifstream, line, '\n')) {
		line = tnstring::trim(line);
		if (line[0] != '#' && line != "") {
				
			numscan = sscanf(line.c_str(), "%s", contact_matrix_filename);
			
			ostr.str("");
			ostr << "+---------------------------------------------------";
			for (unsigned int i = 0; contact_matrix_filename[i] != '\0'; ++i) ostr << "-";
			ostr << "---------------------------------+" << std::endl;
			std::cerr << ostr.str();
			std::cerr << "|                                 Processing File: " << contact_matrix_filename << "                                  |" << std::endl;
			std::cerr << ostr.str();
				
			//Loading contact matrix file to memory
			contact_matrix_ifstream.open(contact_matrix_filename, std::ifstream::in);
			if (contact_matrix_ifstream.fail() || !contact_matrix_ifstream.is_open()) {
				std::cerr << "WARNING: Cannot open file \"" << contact_matrix_ifstream << "\"." << std::endl;
			} else {

				mtx = tnpot::upload_rna_contact_matrix(contact_matrix_ifstream);
				contact_matrix_ifstream.close();
				
				if (n == 0) {
					dim = mtx.get_dimensions();
					sum = new matrix(4, dim[0], dim[1], dim[2], dim[3]);
					pot = new matrix(4, dim[0], dim[1], dim[2], dim[3]);
					n++;
				}
			
				tnpot::sum_rna_contact_matrices(sum, &mtx);
			}
		}
	}
	
	contact_matrix_list_ifstream.close();
	
	//tnpot::print_rna_contact_matrix(*sum, std::cout);
	
	if (t_contact_matrix_file) {
		out.open(t_contact_matrix_filename, std::iostream::out);
		tnpot::print_rna_contact_matrix(*sum, out);
		out.close();
	}

	tnpot::get_rna_potential_matrix(pot, sum);

	if (out_file) {
		out.open(out_filename, std::iostream::out);
		tnpot::print_rna_potential_matrix(*pot, out);
		out.close();
	} else {
		tnpot::print_rna_potential_matrix(*pot, std::cout);
	}
	
	if (print_asgl_plot) tnpot::asgl_potential_plot(pot);

    return 0;
}


void usage(char *prgn){
    std::cout <<
		"Usage: " << prgn << " [-g] [-f <filename>] [-o <filename>] -i <filename>\n"
		"\tParameters:\n"
		"\t-i <filename>        File containing a list of RNA contact matrix paths.\n\n"
		"\tOptons:\n"
		"\t-g                   Print ASGL plots\n"
		"\t-f <filename>        Print Total Contact Matrix\n"
		"\t-o <filename>        File for writing output\n\n"
		"\tRNA POTENTIAL\n"
		"\tThis programme calculates the potential matrix according to a ist of RNA contact matrices\n\n"
		"\tTomas Norambuena A. <tanoramb@puc.cl>\n";
	
	return;
}


