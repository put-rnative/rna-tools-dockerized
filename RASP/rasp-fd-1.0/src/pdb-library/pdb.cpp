#include <algorithm>
#include <cassert>
#include <ctime>
#include <cctype>
#include "../tool-library/string.h"
#include "pdb.h"
#include "format.h"

namespace tnpdb {
	
	//Constructor
	//////////////////////////////////////////////////////////////////////////////////////////////
	pdb::pdb(){
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	pdb::pdb(std::istream &in, bool errors) {
		initialize_hvatm_res_data();
		load(in, errors);
		set_terminal_atoms();
	}

	//Retrieve data
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<atom *> pdb::get_all_atoms(residue::label la) {
		std::vector<atom *>
			atoms;
			
		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues,
			n_atoms;
			
		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					if (_models[i]._chains[j]._residues[k]._label == la) {
						n_atoms = _models[i]._chains[j]._residues[k].get_number_of_atoms();
						for (unsigned int m = 0; m < n_atoms; ++m) {
							atoms.push_back(&_models[i]._chains[j]._residues[k]._atoms[m]);
						}
					}
				}
			}
		}       
	
		return atoms;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<atom *> pdb::get_all_atoms() {
		std::vector<atom *>
			atoms;
			
		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues,
			n_atoms;
			
		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					//if (_models[i]._chains[j]._residues[k]._label == la) {
						n_atoms = _models[i]._chains[j]._residues[k].get_number_of_atoms();
						for (unsigned int m = 0; m < n_atoms; ++m) {
							atoms.push_back(&_models[i]._chains[j]._residues[k]._atoms[m]);
						}
					//}
				}
			}
		}       
	
		return atoms;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<residue *> pdb::get_all_residues(residue::label la) {
		std::vector<residue *>
			residues;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					//std::cerr << residue_name_to_string(_models[i]._chains[j]._residues[k]._name) << " " << _models[i]._chains[j]._residues[k]._label << " " << la << std::endl;
					if (_models[i]._chains[j]._residues[k]._label == la) {
						residues.push_back(&_models[i]._chains[j]._residues[k]);
						//std::cerr << residue_name_to_string(residues.back()->get_name()) << std::endl;
					}
				}
			}
		}
		
		return residues;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<residue *> pdb::get_all_residues() {
		std::vector<residue *>
			residues;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					//if (_models[i]._chains[j]._residues[k]._label == la) {
						residues.push_back(&_models[i]._chains[j]._residues[k]);
					//}
				}
			}
		}
		
		return residues;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<char> pdb::get_all_residue_chain_ids(residue::label la) {
		std::vector<char>
			char_id;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					if (_models[i]._chains[j]._residues[k]._label == la) {
						char_id.push_back(_models[i]._chains[j].get_chain_id());
					}
				}
			}
		}
		
		return char_id;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<char> pdb::get_all_residue_chain_ids() {
		std::vector<char>
			char_id;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					//if (_models[i]._chains[j]._residues[k]._label == la) {
						char_id.push_back(_models[i]._chains[j].get_chain_id());
					//}
				}
			}
		}
		
		return char_id;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	/*std::vector<int> pdb::get_all_residue_indexes(residue::label la) {
		std::vector<int>
			indexes;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					if (_models[i]._chains[j]._residues[k]._label == la) {
						indexes.push_back(int(_models[i]._chains[j]._residues[k].get_index()));
					}
				}
			}
		}
		
		return indexes;
	}*/
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	/*std::vector<int> pdb::get_all_residue_indexes() {
		std::vector<int>
			indexes;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				n_residues = _models[i]._chains[j].get_number_of_residues();
				for (unsigned int k = 0; k < n_residues; ++k) {
					//if (_models[i]._chains[j]._residues[k]._label == la) {
						indexes.push_back(int(_models[i]._chains[j]._residues[k].get_index()));
					//}
				}
			}
		}
		
		return indexes;
	}*/

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<std::string> pdb::get_chain_ids(residue::label la) const {
		std::vector<std::string>
			chains;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				if (_models[i]._chains[j]._residues[0]._label == la) {
					chains.push_back(std::string(1,toupper(_models[i]._chains[j].get_chain_id())));
				}	
			}
		}

		return chains;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<std::string> pdb::get_chain_ids() const {
		std::vector<std::string>
			chains;

		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				chains.push_back(std::string(1,toupper(_models[i]._chains[j].get_chain_id())));	
			}
		}

		return chains;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::string pdb::get_sequence(char chain) const {
		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		std::string
			seq = "";

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				if (_models[i]._chains[j].get_chain_id() == chain && tnstring::trim(std::string(1, _models[i]._chains[j].get_chain_id())) != "") {
					n_residues = _models[i]._chains[j].get_number_of_residues();
					for (unsigned int k = 0; k < n_residues; ++k) {
						seq = seq + residue_name_to_1_letter_string(_models[i]._chains[j]._residues[k].get_name());
					}
					break;
				}
			}
		}

		return seq;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<int> pdb::get_index_sequence(char chain) const {
		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		std::vector<int>
			seq;

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				if (_models[i]._chains[j].get_chain_id() == chain && tnstring::trim(std::string(1, _models[i]._chains[j].get_chain_id())) != "") {
					n_residues = _models[i]._chains[j].get_number_of_residues();
					for (unsigned int k = 0; k < n_residues; ++k) {
						seq.push_back(int(_models[i]._chains[j]._residues[k].get_index()));
					}
					break;
				}
			}
		}

		return seq;
	}
	
	//Add/Set/Modify data
	//////////////////////////////////////////////////////////////////////////////////////////////
	void pdb::add_model(const model &mod) {
		_models.push_back(mod);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	index<residue> pdb::renumber_residues(char chain, index<residue> ini) {
		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues;

		int
			_index = int(ini);

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				if (_models[i]._chains[j].get_chain_id() == chain) {
					n_residues = _models[i]._chains[j].get_number_of_residues();
					for (unsigned int k = 0; k < n_residues; ++k) {
						_models[i]._chains[j]._residues[k].set_index(index<residue>(_index));
						_index++;
					}
					break;
				}
			}
		}

		return index<residue>(_index);

	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	index<atom> pdb::renumber_atoms(char chain, index<atom> ini) {
		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains,
			n_residues,
			n_atoms;

		int
			_index = int(ini);

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned j = 0; j < n_chains; ++j) {
				if (_models[i]._chains[j].get_chain_id() == chain) {
					n_residues = _models[i]._chains[j].get_number_of_residues();
					for (unsigned int k = 0; k < n_residues; ++k) {
						n_atoms = _models[i]._chains[j]._residues[k].get_number_of_atoms();
						for (unsigned int m = 0; m < n_atoms; ++m) {
							_models[i]._chains[j]._residues[k]._atoms[m].set_index(index<atom>(_index));
							_index++;
						}
					}
					break;
				}
			}
		}

		return index<atom>(_index + 1);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	int pdb::rename_chain(char or_chain, char re_chain) {
	
		if (or_chain == re_chain) return -1;

		std::vector<chain>::iterator
			it_chain,
			it_chain_rm,
			it_chain_ex;

		bool
			is_set_it_chain_ex = false;

		unsigned int
			//n_models = _models.size(),
			n_models = 1;

		for (unsigned int i = 0; i < n_models; ++i) {
			for (it_chain = _models[i].chains_begin(); it_chain != _models[i].chains_end(); ++it_chain) {
				if (it_chain->get_chain_id() == re_chain) {
					it_chain_ex = it_chain;
					is_set_it_chain_ex = true;
					break;
				}
			}
			for (it_chain = _models[i].chains_begin(); it_chain != _models[i].chains_end(); ++it_chain) {
				if (it_chain->get_chain_id() == or_chain) {
					it_chain->set_chain_id(re_chain);
					it_chain_rm = it_chain;
					if (is_set_it_chain_ex) {
						if (tnstring::compare_nocase(std::string(1, or_chain), std::string(1, re_chain))) {
							it_chain_ex->_residues.insert(it_chain_ex->residues_begin(), it_chain_rm->residues_begin(), it_chain_rm->residues_end());
						} else {
							it_chain_ex->_residues.insert(it_chain_ex->residues_end(), it_chain_rm->residues_begin(), it_chain_rm->residues_end());
						}
						_models[i]._chains.erase(it_chain_rm);
					}
					break;
				}
			}
		}

		if (is_set_it_chain_ex) return 1;
		else return 0;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void pdb::sort_chains(bool sort_by_label) {
	
		unsigned int
			//n_models = _models.size(),
			n_models = 1,
			n_chains;

		std::vector<chain>
			chains_aa,
			chains_nu,
			chains_xx;
	
		if (sort_by_label) {

			for (unsigned int i = 0; i < n_models; ++i) {
				n_chains = _models[i].get_number_of_chains();
				for (unsigned j = 0; j < n_chains; ++j) {
					if (_models[i]._chains[j]._residues[0]._label == residue::AA) {
						chains_aa.push_back(_models[i]._chains[j]);
					} else if (_models[i]._chains[j]._residues[0]._label == residue::NU) {
						chains_nu.push_back(_models[i]._chains[j]);
					} else {
						chains_xx.push_back(_models[i]._chains[j]);
					}
				}

				_models[i]._chains.clear();

				if (chains_xx.size() > 0) {
					std::sort(chains_xx.begin(), chains_xx.end(), chain_id_less_than());
					_models[i]._chains.insert(_models[i]._chains.begin(), chains_xx.begin(), chains_xx.end());
					chains_xx.clear();
				}
				if (chains_nu.size() > 0) {
					std::sort(chains_nu.begin(), chains_nu.end(), chain_id_less_than());
					_models[i]._chains.insert(_models[i]._chains.begin(), chains_nu.begin(), chains_nu.end());
					chains_nu.clear();
				}
				if (chains_aa.size() > 0) {
					std::sort(chains_aa.begin(), chains_aa.end(), chain_id_less_than());
					_models[i]._chains.insert(_models[i]._chains.begin(), chains_aa.begin(), chains_aa.end());
					chains_aa.clear();
				}
			}
		} else {
			for (unsigned int i = 0; i < n_models; ++i) {
				n_chains = _models[i].get_number_of_chains();
				for (unsigned j = 0; j < n_chains; ++j) {
					if (_models[i]._chains[j]._residues[0]._label == residue::AA || _models[i]._chains[j]._residues[0]._label == residue::NU) {
						chains_aa.push_back(_models[i]._chains[j]);
					} else {
						chains_xx.push_back(_models[i]._chains[j]);
					}
				}

				_models[i]._chains.clear();

				if (chains_xx.size() > 0) {
					std::sort(chains_xx.begin(), chains_xx.end(), chain_id_less_than());
					_models[i]._chains.insert(_models[i]._chains.begin(), chains_xx.begin(), chains_xx.end());
					chains_xx.clear();
				}
				if (chains_aa.size() > 0) {
					std::sort(chains_aa.begin(), chains_aa.end(), chain_id_less_than());
					_models[i]._chains.insert(_models[i]._chains.begin(), chains_aa.begin(), chains_aa.end());
					chains_aa.clear();
				}
			}
		}
		
	}

	//Output
	//////////////////////////////////////////////////////////////////////////////////////////////
	void pdb::write_pdb(std::ostream &out, bool print_header) const {
		unsigned int
			n_hdr = _header.size(),
			n_mod = _models.size();
		
		if (print_header) {
			for (unsigned int i = 0; i < n_hdr; ++i) {
				out << _header[i] << std::endl;
			}
		} else {
			time_t rawtime;
			struct tm *timeinfo;
			time(&rawtime);
			timeinfo = localtime (&rawtime);
			out << "HEADER    PDB FORMAT FILE          " << asctime(timeinfo);
			out << "REMARK   1 PDB FILE MODIFIED" << std::endl;
		}
		
		for (unsigned int i = 0; i < n_mod; ++i) {
			_models[i].write_pdb(out);
		}
		out << "END   \n";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void pdb::write_pdb(std::ostream &out, residue::label la, bool print_header) const {
		unsigned int
			n_hdr = _header.size(),
			n_mod = _models.size();
	
		if (print_header) {
			for (unsigned int i = 0; i < n_hdr; ++i) {
				out << _header[i] << std::endl;
			}
		} else {
			time_t rawtime;
			struct tm *timeinfo;
			time(&rawtime);
			timeinfo = localtime (&rawtime);
			out << "HEADER    PDB FORMAT FILE          " << asctime(timeinfo);
			out << "REMARK   1 PDB FILE MODIFIED" << std::endl;
		}

		for (unsigned int i = 0; i < n_mod; ++i) {
			_models[i].write_pdb(out,la);
		}
		out << "END   \n";
	}


	//Model Iterators
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<model>::const_iterator pdb::models_begin() const {
		return _models.begin();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<model>::const_iterator pdb::models_end() const {
		return _models.end();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<model>::iterator pdb::models_begin() {
		return _models.begin();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<model>::iterator pdb::models_end() {
		return _models.end();
	}


	//Protected
	//////////////////////////////////////////////////////////////////////////////////////////////
	void pdb::load(std::istream &in, bool errors) {
		std::string
			line;

		record
			rec;

		unsigned int
			mnum = 0,
			n_mod = _models.size();
			
		char
			buf[81];

		index<model>
			imdummy(0);
			
		while (getline(in, line, '\n')) {
			rec = record_line(line.c_str());
			if (rec == R_HEADER || rec == R_DBREF || rec == R_SEQRES) {
				_header.push_back(line.c_str());
			} else if (rec == R_MODEL) {
				sscanf(line.c_str(), "%s %d", buf, &mnum);
				
				index<model>
					im(mnum);

				add_model(model(im));
			} else if (rec == R_HETATM || rec == R_ATOM || rec == R_TER || rec == R_ENDMDL) {
				if (_models.empty()) {
					add_model(model(imdummy));
				}
				_models.back().process_record_line(line.c_str());
			} else if (rec == R_MASTER) {
			} else if (rec == R_END) {
			}
		}
		
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void pdb::set_terminal_atoms() {
		unsigned int
			n_models = _models.size(),
			n_chains,
			n;
			
		atom
			*a;
			
		residue::label
			rl_ini,
			rl_end,
			rl;
			
		char
			chn;
			
		bool
			rl_ini_found,
			rl_end_found;
		
		unsigned int
			x = 0, //residue label ini index
			y = 0; //residue label end index

		for (unsigned int i = 0; i < n_models; ++i) {
			n_chains = _models[i].get_number_of_chains();
			for (unsigned int j = 0; j < n_chains; ++j) {
				n = _models[i]._chains[j]._residues.size(); // number of residues
				chn = _models[i]._chains[j].get_chain_id();
				
				if (n > 0) {
					
					rl_ini_found = rl_end_found = false;
					x = y = 0;
					
					for (unsigned int k = 0; k < n; k++) {
						rl = _models[i]._chains[j]._residues[k].get_label();
						
						if (!rl_ini_found && (rl == residue::AA || rl == residue::NU)) {
							rl_ini_found = true;
							x = k;
						}
						
						if (rl_ini_found && !rl_end_found && rl != residue::AA && rl != residue::NU ) {
							rl_end_found = true;
							y = k - 1;
							break;
						}
						
						if (k == n - 1) y = k;
					}
							
					rl_ini = _models[i]._chains[j]._residues[x].get_label();
					rl_end = _models[i]._chains[j]._residues[y].get_label();
					
					//std::cerr << "rINI " << rl_ini << " rEND " << rl_end << std::endl; 
					
					if (rl_ini == residue::AA) {
						a = &_models[i]._chains[j]._residues[x].get_atom_to_modify(atom::AA_N);
						if (a->get_name() != atom::N_INV) a->set_ter(true);
						else std::cerr << "WARNING: Terminal atom " << atom_name_to_string(atom::AA_N, rl_ini) << " in chain " << chn << " not found." << std::endl;
						
						//std::cerr << "INI "<< residue_name_to_string(_models[i]._chains[j]._residues[x].get_name()) << " " << _models[i]._chains[j]._residues[x].get_index()<< std::endl;
						
						//std::cerr << "WARNING: Found terminal atom: ";
						//std::cerr << atom_name_to_string(_models[i]._chains[j]._residues[x].get_atom_to_modify(atom::AA_N).get_name(), _models[i]._chains[j]._residues[x].get_label()) << std::endl;						
					
					} else if (rl_ini == residue::NU) {
						
						a = &_models[i]._chains[j]._residues[x].get_atom_to_modify(atom::NU_OP3);
						if (a->get_name() != atom::N_INV) a->set_ter(true);
						else std::cerr << "WARNING: Terminal atom " << atom_name_to_string(atom::NU_OP3, rl_ini) << " in chain " << chn << " not found." << std::endl;
						
						//_models[i]._chains[j]._residues[0].get_atom_to_modify(atom::NU_O5p).set_ter(true);
						//std::cerr << "WARNING: Found terminal atom: ";
						//std::cerr << atom_name_to_string(_models[i]._chains[j]._residues[x].get_atom_to_modify(atom::NU_O5p).get_name(), _models[i]._chains[j]._residues[0].get_label()) << std::endl;
					
					}
					
					if (rl_end == residue::AA) {
						
						a = &_models[i]._chains[j]._residues[y].get_atom_to_modify(atom::AA_OXT);
						if (a->get_name() != atom::N_INV) a->set_ter(true);
						else std::cerr << "WARNING: Terminal atom " << atom_name_to_string(atom::AA_OXT, rl_end) << " in chain " << chn << " not found." << std::endl;
						
						//std::cerr << "END "<< residue_name_to_string(_models[i]._chains[j]._residues[y].get_name()) << " " << _models[i]._chains[j]._residues[y].get_index()<< std::endl;
						
						//_models[i]._chains[j]._residues[y].get_atom_to_modify(atom::AA_OXT).set_ter(true);
						//std::cerr << "WARNING: Found terminal atom: ";
						//std::cerr << atom_name_to_string(_models[i]._chains[j]._residues[y].get_atom_to_modify(atom::AA_OXT).get_name(), _models[i]._chains[j]._residues[y].get_label()) << std::endl;
					
					} else if (rl_end == residue::NU) {
						
						a = &_models[i]._chains[j]._residues[y].get_atom_to_modify(atom::NU_O3p);
						if (a->get_name() != atom::N_INV) a->set_ter(true);
						else std::cerr << "WARNING: Terminal atom " << atom_name_to_string(atom::NU_O3p, rl_end) << " in chain " << chn << " not found." << std::endl;
						
						//_models[i]._chains[j]._residues[y].get_atom_to_modify(atom::NU_O3p).set_ter(true);
						//std::cerr << "WARNING: Found terminal atom: ";
						//std::cerr << atom_name_to_string(_models[i]._chains[j]._residues[y].get_atom_to_modify(atom::NU_O3p).get_name(), _models[i]._chains[j]._residues[n - 1].get_label()) << std::endl;
					
					}
				}
			}
		}
	}


};
