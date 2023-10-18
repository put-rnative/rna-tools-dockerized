#include "chain.h"
#include "format.h"

namespace tnpdb {

	static residue
		dummy_residue;

	static atom
		dummy_atom;

	// Constructor
	//////////////////////////////////////////////////////////////////////////////////////////////
	chain::chain(){
		_chain = ' ';
	}


	//Retrieve data
	//////////////////////////////////////////////////////////////////////////////////////////////
	char chain::get_chain_id() const {
		return _chain;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int chain::get_number_of_residues() const {
		return _residues.size();
	}
	
	//Set data
	//////////////////////////////////////////////////////////////////////////////////////////////
	void chain::set_chain_id(char c) {
		_chain = c;
	}
	
	//Output
	//////////////////////////////////////////////////////////////////////////////////////////////
	void chain::write_pdb(std::ostream &out) const {
		char
			line[81];

		unsigned int
			n_res = _residues.size();
		
		for (unsigned int i = 0; i < n_res; i++) {
			const residue &res = _residues[i];
			res.write_pdb(_chain, out);
		}

		if (!_residues.empty()) {
			sprintf(line, ter_line_oformat,
				int(_residues.back().get_last_atom_index()) + 1,
				residue_name_to_string(_residues.back().get_name()).c_str(),
				get_chain_id(),
				int(_residues.back().get_index()),
				' '
			);
			out << line << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void chain::write_pdb(std::ostream &out, residue::label la) const {
		char
			line[81];

		unsigned int
			n_res = _residues.size();

		bool
			print_ter = false;

		for (unsigned int i = 0; i < n_res; i++) {
			if (_residues[i]._label == la) { 
				_residues[i].write_pdb(_chain, out);
				print_ter = true;
			}
		}

		if (!_residues.empty() && print_ter) {
			sprintf(line, ter_line_oformat,
				int(_residues.back().get_last_atom_index()) + 1,
				residue_name_to_string(_residues.back().get_name()).c_str(),
				get_chain_id(),
				int(_residues.back().get_index()),
				' '
			);
			out << line << std::endl;
		}
	}

	//Residue Iterators
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<residue>::const_iterator chain::residues_begin() const {
		return _residues.begin();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<residue>::const_iterator chain::residues_end() const {
		return _residues.end();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<residue>::iterator chain::residues_begin() {
		return _residues.begin();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<residue>::iterator chain::residues_end() {
		return _residues.end();
	}

	// Protected
	//////////////////////////////////////////////////////////////////////////////////////////////
	void chain::process_record_line(const char *line, char chn) {
		record
			rec = record_line(line);
		
		if (rec == R_ATOM || rec == R_HETATM) {
			// the read values are not zero padded so we must fill the buffers for strings with 0s.
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
				_chn,
				i_code,
				seg_id[5] = {'\0'},
				element[3] = {'\0'},
				charge[3] = {'\0'},
				_rec;

			//bool
			//	first_res = false,
			//	ter = false;

			std::string
				format;

			if (rec == R_ATOM) {
				format = atom_line_iformat;
				_rec = 'A';
			} else {
				format = hetatm_line_iformat;
				_rec = 'H';
			}
			
			numscan = sscanf(line, format.c_str(),
				&atmindex,
				atmname,
				&alt,
				resname,
				&_chn,
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
			
			/*if (_chain == ' ')*/ _chain = _chn;


			//std::cerr << "CHAIN: " << _chain << std::endl;

			if (/*_chain != ' ' &&*/ chn != _chain) {
				std::cerr << "Confusion over chain labels. Expected " << chn << " got " << _chn << " on line:\n" << line << std::endl;
				return;
			}
			
			if (resindex < 0) {
				std::cerr << "Got negative residue index on line:\n" << line << std::endl;
			//	return;
			}
			
			index<residue>
				_resindex(resindex);

		
			residue::name
				_resname = residue_string_to_name(resname);


			residue::label
				_reslabel = residue_name_to_label(_resname);

			//if (_resname == residue::ZN) std::cerr << "____>" << _resname << " " << _reslabel << std::endl;

			if (_resname == residue::N_INV || _reslabel == residue::T_INV) {
				std::cerr << "Index: " << _resindex << std::endl;
			}

			atom::name
				_atmname = atom_string_to_name(atmname, _reslabel);

			if (_atmname == atom::N_INV && _resname != residue::N_INV && _reslabel != residue::T_INV) {
				std::cerr << "Index: " << atmindex << std::endl;
			}


				
			if (_resname != residue::N_INV && _atmname != atom::N_INV) {
				if (_residues.empty() || _residues.back().get_index() != _resindex) {
					//if(_residues.empty()) first_res = true;
					_residues.push_back(residue(_resname));
					_residues.back().set_index(_resindex);
					_residues.back().set_label(_reslabel);
					_residues.back().set_chain(_chain);
				}
				
				index<atom>
					_atmindex(atmindex);
				
				atom
					atm;
				
				atm.set_coords(point(x,y,z));
				atm.set_index(_atmindex);
				
				if (numscan > 10) {
					atm.set_occupancy(occupancy);
				}
				
				if (numscan > 11) {
					atm.set_temp_factor(temp_factor);
				}
				
				atm.set_segment_id(seg_id);
				atm.set_element(element);
				atm.set_charge(charge);
				atm.set_name(_atmname);
				atm.set_rec(_rec);
				atm.set_alt_loc(alt);
				//if ((_atmname == AA_N && first_res) || _atmname == AA_OXT || (_atmname == NU_O5p &&  
				atm.set_ter(false);
				
				if (int(_residues.back().get_index()) != int(resindex)) {
					std::cerr << "Confusion over residue numbers. Expected " << _residues.back().get_index() << " got " << resindex << " on line:\n" << line << std::endl;
					return;
				}
				
				if (_residues.back().get_name() != residue_string_to_name(resname)) {
					std::cerr << "Confusion over residue types. Expected " << residue_name_to_string(_residues.back().get_name()) << " got " << resname << " on line:\n" << line << std::endl;
					return;
				}
				
				_residues.back().add_atom(atm);

			
			} else {
			}
		} else {
			assert(0);
		}
	}

};

