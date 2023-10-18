#include <vector>
#include <string>
#include "chain.h"

namespace tnpdb {

	class model {
		friend class pdb;
		public:
			// Constructor
			model();
			model(index<model>);

			// Retrieve data
			unsigned int get_number_of_chains() const;
			chain &get_chain(unsigned int);
			const chain &get_chain(unsigned int) const;
			index<model> get_index() const;

			// Set data
			void add_chain(const chain &);
			void set_index(index<model>);

			//Output
			void write_pdb(std::ostream &) const;
			void write_pdb(std::ostream &, residue::label, bool a = false) const;

			//Chain Iterators
			std::vector<chain>::const_iterator chains_begin() const;
			std::vector<chain>::const_iterator chains_end() const;
			std::vector<chain>::iterator chains_begin();
			std::vector<chain>::iterator chains_end();
		
		protected:
    		void process_record_line(const char *);

		private:
			std::vector<chain> _chains;
			std::vector<std::string> _extra;
			index<model> _index;
	};
};

