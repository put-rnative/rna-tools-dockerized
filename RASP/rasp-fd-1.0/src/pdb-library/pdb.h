#include <iostream>
#include <vector>
#include <string>
#include "model.h"

namespace tnpdb {
	
	class pdb {
		public:
			//Constructor
			pdb();
			pdb(std::istream &, bool errors = false);

			//Retrieve data
			std::vector<atom *> get_all_atoms(residue::label);
			std::vector<atom *> get_all_atoms();
			std::vector<residue *> get_all_residues(residue::label);
			std::vector<residue *> get_all_residues();
			std::vector<char> get_all_residue_chain_ids(residue::label);
			std::vector<char> get_all_residue_chain_ids();
			//std::vector<int> get_all_residue_indexes(residue::label);
			//std::vector<int> get_all_residue_indexes();
			std::vector<std::string> get_chain_ids(residue::label) const;
			std::vector<std::string> get_chain_ids() const;
			std::string get_sequence(char) const;
			std::vector<int> get_index_sequence(char) const;

			//Add/Set/Modify data
			void add_model(const model &);
			index<residue> renumber_residues(char, index<residue>);
			index<atom> renumber_atoms(char, index<atom>);
			int rename_chain(char, char);
			void sort_chains(bool a = true);

			//Output
			void write_pdb(std::ostream &, bool a = true) const;
			void write_pdb(std::ostream &, residue::label, bool a = false) const;

			//Model Iterators
			std::vector<model>::const_iterator models_begin() const;
			std::vector<model>::const_iterator models_end() const;
			std::vector<model>::iterator models_begin();
			std::vector<model>::iterator models_end();

		protected:
			void load(std::istream &, bool errors = false);
			void set_terminal_atoms();

		private:
			std::vector<model> _models;
			std::vector<std::string> _header;
	};
};

