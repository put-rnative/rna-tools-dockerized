#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include "../pdb-library/format.h"
#include "../pdb-library/pdb.h"
#include "../tool-library/matrix.h"
#include "../tool-library/string.h"
#include "../tool-library/gradient.h"

#define SIGMA 0.02
#define RT 0.582
#define FLAG_UP 10000000.0
#define FLAG_DOWN -10000000.0

namespace tnpot {

	struct _atom_type {
		tnpdb::residue::label _reslabel;
		tnpdb::residue::name _resname;
		tnpdb::atom::name _atmname;
		bool _terminal;
		bool _centroid;
		unsigned int _ntype;
		std::vector<unsigned int> _type;
	};

	struct _distance_class {
		float _lo_limit;
		float _up_limit;
		unsigned int _class;
	};

	struct _centroid {
		tnpdb::residue::label _reslabel;
		tnpdb::residue::name _resname;
		tnpdb::point _coords;
		unsigned int _n_atm;
		unsigned int _type;
	};

	struct _point_type {
		tnpdb::point _coords;
		unsigned int _type;
		int _resindex;
		int _original_resindex;
		int _original_atmindex;
		std::string _atmname;
		std::string _resname;
		std::string _element;
		char _rec;
		char _alt;
		char _chn;
		char _original_chn;
		bool _ter;
		float _occupancy;
		float _value;
		float _value2;
		float _value3;
	};

	typedef _atom_type atom_type;
	typedef _distance_class distance_class;
	typedef _centroid centroid;
	typedef _point_type point_type;

	struct _interaction {
		point_type *_i;
		point_type *_j;
		float _distance;
		unsigned int _distance_class;
		unsigned int _value;
	};

	typedef _interaction interaction;

	struct _rename_atom {
		tnpdb::residue::label _reslabel;
		tnpdb::residue::name _resname;
		tnpdb::atom::name _atmname_1;
		tnpdb::atom::name _atmname_2;
	};

	typedef _rename_atom rename_atom;
	
	struct _pick {
		int _resindex;
		std::string _resname;
		char _chn;
		float _value;
	};

	typedef _pick pick;

	extern std::vector<atom_type> ATOM_TYPE_DATA;
	extern std::vector<distance_class> DISTANCE_CLASS_DATA;
	extern std::vector<rename_atom> RENAME_ATOM_DATA;

	const char
		atom_type_iformat[] = "%s\t%s\t%s\t%u\t%u\t%s",
		distance_class_iformat[] = "%f\t%f\t%u",
		pdb_list_iformat[] = "%s\t%s",
		rename_atom_iformat[] = "%s\t%s\t%s\t%s";

	void
		load_atom_type_file(std::istream &),
		load_distance_class_file(std::istream &),
		get_rna_contact_matrix(matrix *, std::vector<point_type> *),
		get_rna_potential_matrix(matrix *, matrix *),
		sum_rna_contact_matrices(matrix *, matrix *),
		print_rna_contact_matrix(matrix, std::ostream &),
		print_rna_potential_matrix(matrix, std::ostream &),
		get_nrg_n_contacts(double *, int *, tnpdb::pdb, matrix *),
		get_nrg_mean_nrg_sd(double *, double *, tnpdb::pdb, matrix *, int),
		asgl_potential_plot(matrix *),
		asgl_potential_plots_by_k(std::vector<matrix> *, std::vector<int>),
		print_rasmol_script(std::vector<pick> *);

	unsigned int
		number_of_labels(),
		number_of_atom_types_per_label(tnpdb::residue::label),
		number_of_distance_classes(),
		number_of_centroids(tnpdb::residue::label);

	std::vector<unsigned int>
		get_atom_types_from_index(int);
	
	int
		get_atom_types_index(tnpdb::atom, tnpdb::residue::name),
		get_rename_atom_index(tnpdb::atom, tnpdb::residue::name),
		get_distance_class(float, float a = 2.0);
	
	float
		max_distance();

	std::vector<centroid>
		get_centroids(tnpdb::pdb, tnpdb::residue::label);

	bool
		index_belongs_to_centroid(int),
		is_effective(const point_type &, const point_type &, const std::vector<point_type> &, int, bool a = false, std::ostream &o = std::cout);

	std::vector<point_type>
		get_point_types(tnpdb::pdb, tnpdb::residue::label, bool a = false),
		select_water(tnpdb::pdb, tnpdb::residue::label, float);

	matrix
		contact_matrix_to_diff_table(std::vector<matrix>, int a = 1),
		upload_rna_contact_matrix(std::istream &),
		upload_rna_potential_matrix(std::istream &);
	
	std::vector<pick> 
		get_nrg_profile(tnpdb::pdb, matrix *, int a = 1);

};

