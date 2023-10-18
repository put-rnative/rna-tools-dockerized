#include <cassert>
#include "model.h"
#include "format.h"

namespace tnpdb {
	
	//Constructor
	//////////////////////////////////////////////////////////////////////////////////////////////
	model::model() {
		index<model>
			i;
		_index = i;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	model::model(index<model> i) {
		_index = i;
	}

	// Retrieve data
	//////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int model::get_number_of_chains() const {
		return _chains.size();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	chain &model::get_chain(unsigned int i) {
		assert(i < _chains.size());
		return _chains[i];
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	const chain &model::get_chain(unsigned int i) const {
		 assert(i < _chains.size());
		 return _chains[i];
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	index<model> model::get_index() const {
		return _index;
	}

	// Set data
	//////////////////////////////////////////////////////////////////////////////////////////////
	void model::add_chain(const chain &c) {
		_chains.push_back(c);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void model::set_index(index<model> i) {
		_index = i;
	}

	//Output
	//////////////////////////////////////////////////////////////////////////////////////////////
	void model::write_pdb(std::ostream &out) const {
		char
			line[81];

		unsigned int
			n_chain = _chains.size(),
			n_extra = _extra.size();

		sprintf(line, "MODEL %8d         ", int(_index));
		out << line << std::endl;
		
		for (unsigned int i = 0; i < n_chain; ++i) {
			_chains[i].write_pdb(out);
		}
		
		for (unsigned int i = 0; i < n_extra; ++i) {
			out << _extra[i] << std::endl;
		}
		
		out << "ENDMDL                       " << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void model::write_pdb(std::ostream &out, residue::label la, bool print_extra) const {
		char
			line[81];

		unsigned int
			n_chain = _chains.size(),
			n_extra = _extra.size();

		sprintf(line, "MODEL %8d         ", int(_index));
		out << line << std::endl;

		for (unsigned int i = 0; i < n_chain; ++i) {
			_chains[i].write_pdb(out,la);
		}

		if (print_extra) {
			for (unsigned int i = 0; i < n_extra; ++i) {
				out << _extra[i] << std::endl;
			}
		}

		out << "ENDMDL                       " << std::endl;
	}

	//Chain Iterators
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<chain>::const_iterator model::chains_begin() const {
		return _chains.begin();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<chain>::const_iterator model::chains_end() const {
		return _chains.end();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<chain>::iterator model::chains_begin() {
		return _chains.begin();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<chain>::iterator model::chains_end() {
		return _chains.end();
	}


	//Protected
	//////////////////////////////////////////////////////////////////////////////////////////////
	void model::process_record_line(const char *line) {
		record
	    	rec = record_line(line);

		if (rec == R_ATOM || rec == R_HETATM) {
			int
				atmindex = -1,
				resindex = -1,
				numscan;
						
			float
				x,
				y,
				z,
				occupancy,
				temp_factor;
				
			char
				atmname[5] = {'\0'},
				alt = '\0',
				resname[4] = {'\0'},
				chn,
				i_code,
				seg_id[5] = {'\0'},
				element[3] = {'\0'},
				charge[3] = {'\0'};

			std::string
				format;

			if (rec == R_ATOM) format = atom_line_iformat;
			else format = hetatm_line_iformat;

			numscan = sscanf(line, format.c_str(),
				&atmindex,
				atmname,
				&alt,
				resname,
				&chn,
				&resindex,
				&i_code,
				&x,
				&y,
				&z,
				&occupancy,
				&temp_factor,
				seg_id,
				element,
				charge
			);
				
			//assert(numscan > 5);

			
			if (_chains.empty() || _chains.back().get_chain_id() != chn) {
				_chains.push_back(chain());
			}
				
			_chains.back().process_record_line(line, chn);
																						
		} else if (rec == R_ENDMDL) {
		}
	}

};

