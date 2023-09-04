#include <algorithm>
#include "../tool-library/string.h"
#include "residue.h"
#include "format.h"

namespace tnpdb {
	

	//////////////////////////////////////////////////////////////////////////////////////////////
	residue::residue(residue::name name) {
		_name = name;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	residue::name residue::get_name() const {
		return _name;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	residue::label residue::get_label() const {
		return _label;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	char residue::get_chain() const {
		return _chain;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int residue::get_number_of_atoms() const {
		return _atoms.size();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	const atom &residue::get_atom(atom::name name) const {
		static atom
			dummy_atom;

		//Residue::Atom_label fal= Residue_data::fix_atom_label(label_, al);
		for (std::vector<atom>::const_iterator it = _atoms.begin(); it != _atoms.end(); ++it) {
			if (it->get_name() == name) return *it;
		}
		return dummy_atom;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	const atom &residue::get_atom_at(int i) const {
		static atom
			dummy_atom;

		if (i < 0 || i >= get_number_of_atoms()) {
			std::cerr << "Bad index " << i << " to retrieve atom from residue " << residue_name_to_string(get_name()) << "." << std::endl;
			return dummy_atom;
		}
		return _atoms[i];
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	atom &residue::get_atom_at(int i) {
		static atom
			dummy_atom;

		if (i < 0 || i >= get_number_of_atoms()) {
			std::cerr << "Bad index " << i << " to retrieve atom from residue " << residue_name_to_string(get_name()) << "." << std::endl;
			return dummy_atom;
		}
		return _atoms[i];
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	atom &residue::get_atom_to_modify(atom::name name) {
		static atom
			dummy_atom;
		
		//Residue::Atom_label fal= Residue_data::fix_atom_label(label_, al);
		for (std::vector<atom>::iterator it = _atoms.begin(); it != _atoms.end(); ++it) {
			if (it->get_name() == name) return *it;
		}
		
		return dummy_atom;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	index<atom> residue::get_last_atom_index() const {
		index<atom> max = _atoms.begin()->get_index();
		for (std::vector<atom>::const_iterator it = _atoms.begin(); it != _atoms.end(); ++it) {
			max = std::max(max, it->get_index());
		}
		return max;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	index<residue> residue::get_index() const {
		return _index;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	void residue::set_index(index<residue> ind) {
		_index = ind;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void residue::set_label(residue::label la) {
		_label = la;
	}

    //////////////////////////////////////////////////////////////////////////////////////////////
	void residue::set_chain(char ch) {
		_chain = ch;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void residue::add_atom(const atom &atm) {
		atom::name
			atmname = atm.get_name(); 
		    //Residue::Atom_label al= Residue_data::fix_atom_label(label_, ial);
		if (!can_have_atom(atmname)) {
			std::cerr << "Trying to add invalid atom " <<  atom_name_to_string(atmname, _label) << " on a residue \"" << residue_name_to_string(_name) << "\"" << std::endl;
			return;
		}
										
		if (atmname == atom::N_INV) {
			return;
		}
		
		_atoms.push_back(atm);
		_atoms.back().set_etype(atom_name_to_etype(atmname, _label));
		
		if (_min_atom_index) {
			_min_atom_index = std::min(_min_atom_index, atm.get_index());
		} else {
			_min_atom_index = atm.get_index();
		}
		
	}

	//Output

	//////////////////////////////////////////////////////////////////////////////////////////////
	void residue::write_pdb(char chn, std::ostream &out) const {
		char
			line[81];


		std::vector<atom>
			atoms(_atoms.begin(), _atoms.end());

		unsigned int
			n_atm = atoms.size();


		if (n_atm > 1) std::sort(atoms.begin(), atoms.end(), atom_index_less_than());

		for (unsigned int i = 0; i < n_atm; i++) {
			atom::name
				name = atoms[i].get_name();

			//std::cerr << "atom coords: " << atoms[i].get_coords() << std::endl;
			//std::cerr << "atom name:" << name << std::endl;

			const atom
				&atm = atoms[i];

			point
				pt = atm.get_coords();
			
			char
				alt = atm.get_alt_loc(),
				i_code = ' ',
				rec = atm.get_rec();

			std::string
				format;
			

			if (rec == 'A') format = atom_line_oformat;
			else format = hetatm_line_oformat;


			sprintf(line, format.c_str(),
				int(atm.get_index()),
				atom_name_to_string(name, _label).c_str(),
				alt,
				residue_name_to_string(get_name()).c_str(),
				chn,
				int(get_index()),
				i_code,
				pt.x(),
				pt.y(),
				pt.z(),
				atm.get_occupancy(),
				atm.get_temp_factor(),
				atm.get_segment_id(),
				atm.get_element(),
				atm.get_charge()
			);
			out << line << std::endl;
		}
	}

	//Protected

	//////////////////////////////////////////////////////////////////////////////////////////////
	index<atom> residue::get_min_atom_index() const {
		return _min_atom_index;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	bool residue::can_have_atom(atom::name atmname) const {
		atom::name
			*list_hv,
			*list_hy,
			*list_hvbb,
			*list_hybb;

		bool
			found = false;

		std::vector<atom::name>
			dummy;


		list_hv = hvatm_res_data[_name];
		//list_hy = hyatm_res_data[_name];

		if(!list_hv || !list_hy) {
			std::cerr << "Residue not initialized yet." << std::endl;
			return false;
		}

		if (_name == ATP) {
			list_hvbb = hvatm_dum;
		} else if (_name == ADE || _name == CYT || _name == GUA || _name == THY || _name == URA) {
			list_hvbb = hvatm_res_data[NUB];
			//list_hybb = hyatm_res_data[NUB];
		} else { 
			list_hvbb = hvatm_res_data[AAB];
			//list_hybb = hyatm_res_data[AAB];
		}

		for (unsigned i = 0; list_hvbb[i] != atom::N_INV && !found; ++i) {
			if(list_hvbb[i] == atmname) found = true;
		}

		for (unsigned i = 0; list_hv[i] != atom::N_INV && !found; ++i) {
			if(list_hv[i] == atmname) found = true;
		}

		/*
		for (unsigned i = 0; list_hybb[i] != atom::N_INV && !found; ++i) {
			if(list_hybb[i] == atmname) found = true;
		}

		for (unsigned i = 0; list_hy[i] != atom::N_INV && !found; ++i) {
			if(list_hy[i] == atmname) found = true;
		}*/

		if(found) return true;
		return false;
	}


	//Conversion
	//////////////////////////////////////////////////////////////////////////////////////////////
	atom::etype atom_name_to_etype(atom::name name, residue::label la) {
		if (la == residue::T_INV) return atom::T_INV;	
		atom_data
			*data = get_atom_name_data(la);

		for (unsigned int i = 0; data[i].name != atom::N_INV; ++i) {
			if (name == data[i].name) {
				return data[i].etype;
			}
		}
		std::cerr << "Unknown atom name: " << name << std::endl;
		return atom::T_INV;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	std::string atom_name_to_string(atom::name name, residue::label la) {
		if (la == residue::T_INV) return "UNKN";

		atom_data
			*data = get_atom_name_data(la);

		for (unsigned int i = 0; data[i].name != atom::N_INV; ++i) {
			if (name == data[i].name) {
				return data[i].s;
			}
		}
		std::cerr << "Unknown atom name: " << name << ", returning UNKN" << std::endl;
		return "UNKN";
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	atom::name atom_string_to_name(const std::string &c, residue::label la) {
		if (la == residue::T_INV) return atom::N_INV;

		atom_data
			*data = get_atom_name_data(la);
			
		std::string
			s = tnstring::trim(c);

		for (unsigned int i = 0; data[i].name != atom::N_INV; i++) {
			//if (la == residue::ME) std::cerr << s << " " << tnstring::trim(data[i].s) << std::endl;
			if (s == tnstring::trim(data[i].s)) {
				return data[i].name;
			}
		}
		std::cerr << "\""<< s << "\" is not a known atom name of label " << la << std::endl;
		
		return atom::N_INV;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	residue::name residue_string_to_name(const std::string &c) {
		std::string s;
		s = tnstring::trim(c);
		for (unsigned int i = 0; res_name_data[i].name != residue::N_INV; i++) {
			//if (res_name_data[i].label == residue::ME) std::cerr << "--->\""<< s << "\" " << tnstring::trim(res_name_data[i].s) << std::endl;
			if (s == tnstring::trim(res_name_data[i].s)) {
				//if (res_name_data[i].label == residue::ME) std::cerr << "--->\""<< s << "\" " << tnstring::trim(res_name_data[i].s) << " " << res_name_data[i].name << std::endl;
				return res_name_data[i].name;
			}
		}
		std::cerr << "\""<< s << "\" is not a known residue" << std::endl;;
		return residue::N_INV;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::string residue_name_to_string(residue::name name) {
		for (unsigned int i = 0; res_name_data[i].name != residue::N_INV; i++) {
			if (name == res_name_data[i].name) {
				return res_name_data[i].s;
			}
		}
		std::cerr << "Unknown residue name: " << name << ", returning UKN." << std::endl;
		return "UKN";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	residue::label residue_name_to_label(residue::name name) {
		for (unsigned int i = 0; res_name_data[i].name != residue::N_INV; i++) {
			if (name == res_name_data[i].name) {
				//std::cerr << "*****>" << name << " " << res_name_data[i].label << std::endl;
				return res_name_data[i].label;
			}
		}
		std::cerr << "Unknown residue name: " << name << ", returning invalid type ";
		return residue::T_INV;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::string residue_name_to_1_letter_string(residue::name name) {
		for (unsigned int i = 0; res_1_letter_name_data[i].name != residue::N_INV; i++) {
			if (name == res_1_letter_name_data[i].name) {
				return res_1_letter_name_data[i].s;
			}
		}
		std::cerr << "Unknown residue name: " << name << ", returning invalid/unknown 1-letter name '-'" << std::endl;
		return "-";
	}



	//Other
	//////////////////////////////////////////////////////////////////////////////////////////////
	atom_data *get_atom_name_data(residue::label la) {
		switch (la) {
			case residue::AA:
				return atom_name_aa_data;
			case residue::NU:
				return atom_name_nu_data;
			case residue::WA:
				return atom_name_wa_data;
			case residue::ME:
				return atom_name_me_data;
			default:
				return NULL;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void initialize_hvatm_res_data() {
		hvatm_res_data[residue::ALA] = hvatm_ala;
		hvatm_res_data[residue::ARG] = hvatm_arg;
		hvatm_res_data[residue::ASN] = hvatm_asn;
		hvatm_res_data[residue::ASP] = hvatm_asp;
		hvatm_res_data[residue::CYS] = hvatm_cys;
		hvatm_res_data[residue::GLN] = hvatm_gln;
		hvatm_res_data[residue::GLU] = hvatm_glu;
		hvatm_res_data[residue::GLY] = hvatm_gly;
		hvatm_res_data[residue::HIS] = hvatm_his;
		hvatm_res_data[residue::ILE] = hvatm_ile;
		hvatm_res_data[residue::LEU] = hvatm_leu;
		hvatm_res_data[residue::LYS] = hvatm_lys;
		hvatm_res_data[residue::MET] = hvatm_met;
		hvatm_res_data[residue::PHE] = hvatm_phe;
		hvatm_res_data[residue::PRO] = hvatm_pro;
		hvatm_res_data[residue::SER] = hvatm_ser;
		hvatm_res_data[residue::THR] = hvatm_thr;
		hvatm_res_data[residue::TRP] = hvatm_trp;
		hvatm_res_data[residue::TYR] = hvatm_tyr;
		hvatm_res_data[residue::VAL] = hvatm_val;
		hvatm_res_data[residue::ADE] = hvatm_ade;
		hvatm_res_data[residue::CYT] = hvatm_cyt;
		hvatm_res_data[residue::GUA] = hvatm_gua;
		hvatm_res_data[residue::THY] = hvatm_thy;
		hvatm_res_data[residue::URA] = hvatm_ura;
		hvatm_res_data[residue::AAB] = hvatm_aab;
		hvatm_res_data[residue::NUB] = hvatm_nub;
		hvatm_res_data[residue::HOH] = hvatm_hoh;
		hvatm_res_data[residue::ATP] = hvatm_atp;
		hvatm_res_data[residue::CEA] = hvatm_cea;
		hvatm_res_data[residue::CEN] = hvatm_cen;
		hvatm_res_data[residue::MG] = hvatm_mg;
		hvatm_res_data[residue::ZN] = hvatm_zn;
		hvatm_res_data[residue::CA] = hvatm_ca;
		hvatm_res_data[residue::NA] = hvatm_na;
		hvatm_res_data[residue::NI] = hvatm_ni;
		hvatm_res_data[residue::CL] = hvatm_cl;
		hvatm_res_data[residue::FE] = hvatm_fe;
		
		return;
	}

};

