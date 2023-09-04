#include <iostream>
#include <vector>
#include <iterator>
#include <cassert>
#include "residue.h"


namespace tnpdb {

	class chain {
		friend class pdb;
		friend class model;
		public:
			//Constructor
			chain();
			//Retrieve data
    		char get_chain_id() const;
			unsigned int get_number_of_residues() const;
			//Set data
			void set_chain_id(char);
			//Output
			void write_pdb(std::ostream &) const; 		    				         // Write PDB ATOM record
			void write_pdb(std::ostream &, residue::label) const;

			//Residue Iterators
			std::vector<residue>::const_iterator residues_begin() const;
			std::vector<residue>::const_iterator residues_end() const;
			std::vector<residue>::iterator residues_begin();
			std::vector<residue>::iterator residues_end();

		protected:
			void process_record_line(const char *, char);
		
		private:
			std::vector<residue> _residues;
			char _chain;
	};

	//Comparison class
	class chain_id_less_than {
		public:
			inline bool operator()(const chain &c1, const chain &c2){
				return (tolower(c1.get_chain_id()) < tolower(c2.get_chain_id()));
			}
	};

};

