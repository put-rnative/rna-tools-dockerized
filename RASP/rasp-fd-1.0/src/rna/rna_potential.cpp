#include <map>
#include <cmath>
#include <cstdlib>
#include <iomanip> 
#include "rna_potential.h"

int randomNum(int a, int b){
    return (int)((double)a + (((double)b - (double)a + 1.0)*(double)rand())/((double)RAND_MAX + 1.0));
}

namespace tnpot {
	
	std::vector<atom_type> ATOM_TYPE_DATA;
	std::vector<distance_class> DISTANCE_CLASS_DATA;
	std::vector<rename_atom> RENAME_ATOM_DATA;

	//////////////////////////////////////////////////////////////////////////////////////////////
	void load_atom_type_file(std::istream &in) {
		std::string
			line,
			labelstr;

		std::vector<std::string>
			s_type;

		atom_type
			*data;

		tnpdb::residue::label
			reslab;

		int
			numscan,
			terminal,
			centroid,
			type;

		char
			label[3],
			resname[4],
			atmname[5],
			types[100];

		while (getline(in, line, '\n')) {
			line = tnstring::trim(line);
			if (line[0] != '#' && line != "") {
				numscan = sscanf(line.c_str(), atom_type_iformat,
					label,
					resname,
					atmname,
					&terminal,
					&centroid,
					types
				);

				labelstr = label;

				if(labelstr == "AA") reslab = tnpdb::residue::AA;
				else if(labelstr == "NU") reslab = tnpdb::residue::NU;
				else if(labelstr == "WA") reslab = tnpdb::residue::WA;
				
				data = new atom_type;
				data->_reslabel = reslab;
				data->_resname = tnpdb::residue_string_to_name(resname);
				data->_atmname = tnpdb::atom_string_to_name(atmname, reslab);
				data->_terminal = (bool)terminal;
				data->_centroid = (bool)centroid;
				s_type = tnstring::explode(",",(std::string)types);
				data->_ntype = (unsigned int)s_type.size();

				for (unsigned int i = 0; i < data->_ntype; ++i) {
					data->_type.push_back((unsigned int)tnstring::string2int(s_type[i]));
				}


				ATOM_TYPE_DATA.push_back(*data);
				delete data;
			}
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void load_distance_class_file(std::istream &in) {
		std::string
			line;

		distance_class
			*data;

		tnpdb::residue::label
			reslab;

		int
			numscan,
			_class;

		float
			lo_limit,
			up_limit;

		while (getline(in, line, '\n')) {
			line = tnstring::trim(line);
			if (line[0] != '#' && line != "") {
				numscan = sscanf(line.c_str(), distance_class_iformat,
					&lo_limit,
					&up_limit,
					&_class
				);

				data = new distance_class;
				data->_lo_limit = lo_limit;
				data->_up_limit = up_limit;
				data->_class = _class;

				DISTANCE_CLASS_DATA.push_back(*data);
				delete data;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int number_of_labels() {
		unsigned int
			size = (unsigned int)ATOM_TYPE_DATA.size(),
			groups = 1;

		if (size == 0) return 0;

		tnpdb::residue::label
			reslab = ATOM_TYPE_DATA[0]._reslabel;

		for (unsigned int i = 0; i < size; ++i) {
			if (reslab != ATOM_TYPE_DATA[i]._reslabel) groups++;
			reslab = ATOM_TYPE_DATA[i]._reslabel;
		}

		return groups;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int number_of_atom_types_per_label(tnpdb::residue::label la) {
		unsigned int
			size = (unsigned int)ATOM_TYPE_DATA.size(),
			max = 0;
		

		for (unsigned int i = 0; i < size; ++i) {
			if (ATOM_TYPE_DATA[i]._reslabel == la) {
				for (unsigned int j = 0; j < ATOM_TYPE_DATA[i]._type.size(); ++j) {
					if (max < ATOM_TYPE_DATA[i]._type[j]) max = ATOM_TYPE_DATA[i]._type[j];
				}
			}		
		}

		return max;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int number_of_distance_classes() {
		return DISTANCE_CLASS_DATA.size();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	float  max_distance() {
		return DISTANCE_CLASS_DATA[DISTANCE_CLASS_DATA.size() - 1]._up_limit;;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned int> get_atom_types_from_index(int i) {
		std::vector<unsigned int> dummy;
		if (i < 0 || i >= ATOM_TYPE_DATA.size()) return (dummy);
		return ATOM_TYPE_DATA[i]._type;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	int get_atom_types_index(tnpdb::atom atm, tnpdb::residue::name nm) {
		unsigned int
			size = (unsigned int)ATOM_TYPE_DATA.size();
							        
		for (unsigned int i = 0; i < size; ++i) {
			if (ATOM_TYPE_DATA[i]._resname == nm && ATOM_TYPE_DATA[i]._atmname == atm.get_name() && ATOM_TYPE_DATA[i]._terminal == atm.get_ter()) return i;
		}
			
		return (-1);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	int get_distance_class(float d, float min) {
		unsigned int
			size = (unsigned int)DISTANCE_CLASS_DATA.size();
							        
		for (unsigned int i = 0; i < size; ++i) {
			if (d >= min && DISTANCE_CLASS_DATA[i]._lo_limit <= d && d < DISTANCE_CLASS_DATA[i]._up_limit) return DISTANCE_CLASS_DATA[i]._class; 
		}

		return (-1);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int number_of_centroids(tnpdb::residue::label la) {
		unsigned int
			size = (unsigned int)ATOM_TYPE_DATA.size(),
			num = 0;
		int
			type = -1;

		for (unsigned int i = 0; i < size; ++i) {
			if (ATOM_TYPE_DATA[i]._reslabel == la && ATOM_TYPE_DATA[i]._centroid) {
				for (unsigned int j = 0; j < ATOM_TYPE_DATA[i]._type.size(); ++j) {
					if (type != ATOM_TYPE_DATA[i]._type[j]) {
						type = ATOM_TYPE_DATA[i]._type[j];
						num++;
					}
				}
			}
		}
		return num;		
	}

	//////////////////////////////////////////////////////////////////////////////////////////////  VEr asunto de ALT LOC con Centroids
	std::vector<centroid> get_centroids(tnpdb::pdb in, tnpdb::residue::label la) {
		std::vector<tnpdb::residue *>
			residues = in.get_all_residues(la);

		std::vector<centroid>
			cen;

		centroid
			*c;

		tnpdb::residue::name
			resname;

		tnpdb::index<tnpdb::residue>
			prev_resindex,
			resindex;

		tnpdb::atom
			*atm;

		unsigned int
			n_res = residues.size(),
			n_atm,
			n_cen;

		int
			type = -1,
			ind;

		std::map<int, centroid>
			map_cen;

		std::map<int,centroid>::iterator
			it;


		for (unsigned int i = 0; i < n_res; ++i) {
			n_atm = residues[i]->get_number_of_atoms();
			resname = residues[i]->get_name();
			resindex = residues[i]->get_index();
			for (unsigned int j = 0; j < n_atm; ++j) {
				ind = get_atom_types_index(residues[i]->get_atom_at(j),resname);
				if(ind >= 0) {
					if (ATOM_TYPE_DATA[ind]._centroid) {
						for (unsigned int k = 0; k < ATOM_TYPE_DATA[ind]._type.size(); ++k) {
							type = (int)ATOM_TYPE_DATA[ind]._type[k];
							it = map_cen.find(type);
							if (prev_resindex != resindex || it == map_cen.end()) {
								if (prev_resindex != resindex && !map_cen.empty()) {
									for (it = map_cen.begin(); it != map_cen.end(); it++) {
										cen.push_back(it->second);
									}
									map_cen.clear();
								}
								prev_resindex = resindex;
								c = new centroid;
								c->_reslabel = la;
								c->_resname = resname;
								c->_coords = residues[i]->get_atom_at(j).get_coords();
								c->_type = ATOM_TYPE_DATA[ind]._type[k];
								c->_n_atm = 1;
								map_cen[type] = *c;
							} else {
								map_cen[type]._coords = map_cen[type]._coords + residues[i]->get_atom_at(j).get_coords();
								map_cen[type]._n_atm = map_cen[type]._n_atm + 1;
							}
						}
					}
				}
			}
		}

		for (it = map_cen.begin(); it != map_cen.end(); it++) {
			cen.push_back(it->second);
		}

		map_cen.clear();

		n_cen = cen.size();

		for (unsigned int i = 0; i < n_cen; ++i) {
			cen[i]._coords = cen[i]._coords / (double)cen[i]._n_atm;
		}

		return cen;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	bool index_belongs_to_centroid(int i) {
		if (i < 0 ) return false;
		return ATOM_TYPE_DATA[i]._centroid;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<point_type> get_point_types(tnpdb::pdb in, tnpdb::residue::label la, bool iwtrl) {
		std::vector<point_type>
			pt,
			pt_alt_loc;

		point_type
			*p;
	

		std::vector<tnpdb::residue *>
			residues = in.get_all_residues(la),
			waters;


		std::vector<centroid>
			centroids;

		unsigned int
			n_residues,
			n_atoms,
			n_centroids;

		int
			ind;

		char
			chn = (la == tnpdb::residue::NU)?'A':'B',
			o_chn;

		tnpdb::residue::name
			resname;

		tnpdb::residue::label
			reslabel;

		tnpdb::index<tnpdb::residue>
			resindex;
			
		bool
			alt_loc_found = false;
			
		float
			max_occ;
			
		int
			ind_max_occ = 0;

		
		n_residues = residues.size();
		n_centroids = number_of_centroids(la);

		if (n_centroids > 0) {
			centroids = get_centroids(in,la);
		}

		// ATOMS
		for (unsigned int i = 0; i < n_residues; ++i) {
			n_atoms = residues[i]->get_number_of_atoms();
			resname = residues[i]->get_name();
			reslabel = residues[i]->get_label();
			resindex = residues[i]->get_index();
			o_chn = residues[i]->get_chain();
			for (unsigned int j = 0; j < n_atoms; ++j) {
				ind = get_atom_types_index(residues[i]->get_atom_at(j),resname);
				if (ind >= 0 && !index_belongs_to_centroid(ind)) {
					p = new point_type;
					p->_coords = residues[i]->get_atom_at(j).get_coords();
					p->_type = get_atom_types_from_index(ind)[0];
					p->_resindex = (int)resindex;
					p->_atmname = tnpdb::atom_name_to_string(residues[i]->get_atom_at(j).get_name(), reslabel);
					p->_resname = tnpdb::residue_name_to_string(resname);
					p->_element = residues[i]->get_atom_at(j).get_element();
					p->_rec = residues[i]->get_atom_at(j).get_rec();
					p->_original_resindex = (int)resindex;
					p->_original_atmindex = (int)residues[i]->get_atom_at(j).get_index();
					p->_chn = chn;
					p->_original_chn = o_chn;
					p->_alt = residues[i]->get_atom_at(j).get_alt_loc();
					p->_ter = residues[i]->get_atom_at(j).get_ter();
					p->_occupancy = residues[i]->get_atom_at(j).get_occupancy();
					p->_value = 0.0;
					p->_value2 = 0.0;
					p->_value3 = 0.0;
					
					if (!pt.empty() && p->_atmname == pt.back()._atmname && p->_resindex == pt.back()._resindex && p->_alt != pt.back()._alt) {
						alt_loc_found = true;
						pt_alt_loc.push_back(*p);
					}
					
					//if (!pt.empty()) std::cerr << "---" << p->_atmname << "---" << pt.back()._atmname << "---" << p->_resindex << "---" << pt.back()._resindex << "---" << alt_loc_found << std::endl;
					
					if (!pt.empty() && alt_loc_found && ((p->_resindex == pt.back()._resindex && p->_atmname != pt.back()._atmname) || p->_resname != pt.back()._resname)) {
						alt_loc_found = false;
						max_occ = pt.back()._occupancy;
						for (unsigned int n = 0; n < pt_alt_loc.size(); n++) {
							if (max_occ < pt_alt_loc[n]._occupancy) {
								max_occ = pt_alt_loc[n]._occupancy;
								ind_max_occ = n;
							}
						}
						if (max_occ != pt.back()._occupancy) {
							pt.pop_back();
							pt.push_back(pt_alt_loc[ind_max_occ]);
							pt_alt_loc.clear();
						}
					}
					
					if (!alt_loc_found) pt.push_back(*p);
				}
			}
		}

		// CENTROIDS
		for (unsigned int i = 0; i < n_centroids; ++i) {
			p = new point_type;
			p->_coords = centroids[i]._coords;
			p->_type = centroids[i]._type;
			p->_resindex = pt.back()._resindex + 1;
			p->_original_resindex = pt.back()._resindex + 1;;
			p->_original_atmindex = pt.back()._original_atmindex + 1;
			p->_atmname = " CEN";
			p->_resname = "CEN";
			p->_element = "O";
			p->_rec = 'H';
			p->_rec = chn;
			p->_alt = ' ';
			p->_ter = false;
			p->_value = 0.0;
			p->_value2 = 0.0;
			p->_value3 = 0.0;
			pt.push_back(*p);
		}

		if (iwtrl) {
			// WATERS
			waters = in.get_all_residues(tnpdb::residue::WA);
			n_residues = waters.size();
			for (unsigned int i = 0; i < n_residues; ++i) {
				n_atoms = waters[i]->get_number_of_atoms();
				resindex = waters[i]->get_index();
				for (unsigned int j = 0; j < n_atoms; ++j) {
					ind = get_atom_types_index(waters[i]->get_atom_at(j), tnpdb::residue::HOH);
					if (ind >= 0 && !index_belongs_to_centroid(ind)) {
						p = new point_type;
						p->_coords = waters[i]->get_atom_at(j).get_coords();
						p->_type = get_atom_types_from_index(ind)[0];
						p->_resindex = pt.back()._resindex + 1;
						p->_atmname = " O  ";
						p->_resname = "HOH";
						p->_element = "O";
						p->_rec = 'H';
						p->_original_resindex = (int)resindex;
						p->_original_atmindex = (int)waters[i]->get_atom_at(j).get_index();
						p->_chn = ' ';
						p->_alt = ' ';
						p->_ter = false;
						p->_value = 0.0;
						p->_value2 = 0.0;
						p->_value3 = 0.0;
						pt.push_back(*p);
					}
				}
			}
		}

		return pt;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void get_rna_contact_matrix(matrix *mtx, std::vector<point_type> *pt) {
		
		std::vector<point_type>::iterator i, j;
		
		float
			d = -1.0;
			
		int
			d_class = -1,
			val = 0,
			k_max = mtx->get_dimensions()[0] - 1,
			k = 0;
		
		
		for (i = pt->begin(); i != pt->end(); i++) {
			for (j = pt->begin(); j != pt->end(); j++) {
				
				d = tnpdb::distance(i->_coords, j->_coords);
				d_class = get_distance_class(d);
				
				if (d_class > 0) {
					
					//k = fabs(j->_resindex - i->_resindex) - 1; //Symmetric
					k = j->_resindex - i->_resindex - 1; //Asymmetric
					
					//std::cerr << "resi: "<<i->_resindex << " typi: "<< i->_type - 1 << " resj: " << j->_resindex << " typj: " << j->_type - 1 << " dist: " << d << " dclass: " << d_class - 1 << std::endl;
					//std::cerr << i->_resindex << " "<< i->_type - 1 << " " << j->_resindex << " " << j->_type - 1 << " " << d << " " << i->_occupancy << " " << j->_occupancy << std::endl;
					
					if (k >= 0 && i->_original_chn != j->_original_chn) {
						k = k_max;
						//	std::cerr << i->_original_chn << " " << j->_original_chn << std::endl;
					}
					
					if (k >= 0 && k < k_max) {
					
						val = (int)mtx->get(k, i->_type - 1, j->_type - 1, d_class - 1);
						mtx->put(k, i->_type - 1, j->_type - 1, d_class - 1, val + 1.0);
					
					} else if (k >= k_max){
						//std::cerr << "resi: "<<i->_resindex << " typi: "<< i->_type - 1 << " resj: " << j->_resindex << " typj: " << j->_type - 1 << " dist: " << d << " dclass: " << d_class - 1 << std::endl;
						
						val = (int)mtx->get(k_max, i->_type - 1, j->_type - 1, d_class - 1);
						mtx->put(k_max, i->_type - 1, j->_type - 1, d_class - 1, val + 1.0);
					
					}
				}
			}
		}
		
		return;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void get_rna_potential_matrix(matrix *mtx, matrix *freq) {
		
		std::vector<int>
			dim = freq->get_dimensions();
		
		matrix
			num_fxx(2, dim[0], dim[3]); //k, distance

		std::vector<float>
			den_fxx(dim[0],0.0);
			
		float
			Fijd,
			Mijv,
			num_fxxv,
			den_fxxv,
			nrg;

		matrix
			Mij(3,dim[0],dim[1], dim[2]);

		for (unsigned int k = 0; k < dim[0]; k++) {
			for (unsigned int i = 0; i < dim[1]; i++) {
				for (unsigned int j = 0; j < dim[2]; j++) {
					for (unsigned int d = 0; d < dim[3]; d++) {
						Fijd = freq->get(k, i, j, d);
						Mij.put(k, i, j, Mij.get(k, i, j) + Fijd);
						num_fxx.put(k, d, num_fxx.get(k, d) + Fijd);
						den_fxx[k] = den_fxx[k] + Fijd;
					}
				}
			}
		}
		
		for (unsigned int k = 0; k < den_fxx.size(); k++) {
			if (den_fxx[k] == 0.0) {
				std::cerr << "WARNING: Empty frequency matrix for k = " << k << std::endl;
			}
		}

		for (unsigned int k = 0; k < dim[0]; k++) {
			for (unsigned int i = 0; i < dim[1]; i++) {
				for (unsigned int j = 0; j < dim[2]; j++) {
					for (unsigned int d = 0; d < dim[3]; d++) {
					
						Fijd = freq->get(k, i, j, d);
						//std::cerr << Fijd << std::endl;
						Mijv = Mij.get(k, i, j);
						//std::cerr << Mijv << std::endl;
						num_fxxv = num_fxx.get(k, d);
						//std::cerr << num_fxxv << std::endl;
						den_fxxv = den_fxx[k];
						//std::cerr << den_fxxv << std::endl;
									
						if (Mijv > 0.0 && num_fxxv > 0.0) {

							nrg = RT*(log(1.0 + Mijv*SIGMA) - log(1.0 + Mijv*SIGMA*((Fijd/Mijv)/(num_fxxv/den_fxxv))));
						
							mtx->put(k, i, j, d, nrg);
							
							//std::cerr << "A " << k << " " << i << " " << j << " " << d << " Fij:"<< Fijd << " mij:"<< Mijv << " div:" << Fijd/Mijv << " num:" << num_fxxv << " den:" << den_fxxv << " div:" << num_fxxv/den_fxxv << " nrg:" << nrg << std::endl;
					
						} else {

							nrg = RT*(log(1.0 + Mijv*SIGMA));
						
							mtx->put(k, i, j, d, nrg);
							
							//std::cerr << "B " << k << " " << i << " " << j << " " << d << " Fij:"<< Fijd << " mij:"<< Mijv << " div:" << Fijd/Mijv << " num:" << num_fxxv << " den:" << den_fxxv << " div:" << num_fxxv/den_fxxv << "nrg:" << nrg << std::endl;
					
						}
					}
				}
			}
		}
		
		return;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void print_rna_contact_matrix(matrix mtx, std::ostream &out) {

		std::vector<int>
			dim = mtx.get_dimensions();
	
		out << "# Contact Matrix" << std::endl;
		out << dim[0] << "\t" << dim[1] << "\t" << dim[2] << "\t" << dim[3] << std::endl;
	    out << "# K\tTYPE\tTYPE\tDIST\tFREQ" << std::endl;

		for (unsigned int k = 0; k < dim[0]; ++k) {
			for (unsigned int i = 0; i < dim[1]; ++i) {
				for (unsigned int j = 0; j < dim[2]; ++j) {
					for (unsigned int d = 0; d < dim[3]; ++d) {					
						out << k << "\t"<< i << "\t" << j << "\t" << d << "\t" << (int)mtx.get(k,i,j,d) << std::endl;
					}
				}
			}
		}
		
		return;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	void print_rna_potential_matrix(matrix mtx, std::ostream &out) {
		
		std::vector<int>
			dim = mtx.get_dimensions();
	
		out << "# Potential Matrix" << std::endl;
		out << dim[0] << "\t" << dim[1] << "\t" << dim[2] << "\t" << dim[3] << std::endl;
	    out << "# K\tTYPE\tTYPE\tDIST\tENERGY" << std::endl;

		for (unsigned int k = 0; k < dim[0]; ++k) {
			for (unsigned int i = 0; i < dim[1]; ++i) {
				for (unsigned int j = 0; j < dim[2]; ++j) {
					for (unsigned int d = 0; d < dim[3]; ++d) {					
						out << k << "\t"<< i << "\t" << j << "\t" << d << "\t" << std::setprecision(4) << std::setiosflags(std::ios::fixed)<< mtx.get(k,i,j,d) << std::endl;
					}
				}
			}
		}
		
		return;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	void sum_rna_contact_matrices(matrix *sum, matrix *mtx) {
		std::vector<int>
			dims = sum->get_dimensions(),
			dimm = mtx->get_dimensions();


		for (unsigned int i = 0; i < 4; i++) {
			if (dims[i] != dimm[i]) {
				std::cerr << "ERROR: Cannot sum RNA contact matrices of different dimensions" << std::endl;
				exit(1);
			}
		}	

		for (unsigned int k = 0; k < dims[0]; ++k) {
			for (unsigned int i = 0; i < dims[1]; ++i) {
				for (unsigned int j = 0; j < dims[2]; ++j) {
					for (unsigned int d = 0; d < dims[3]; ++d) {
						sum->put(k, i, j, d, sum->get(k, i, j, d) + mtx->get(k, i, j, d));
					}
				}
			}
		}
		
		return;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	matrix upload_rna_contact_matrix(std::istream &in) {
		matrix
			*mtx;

		bool
			firstline = true;

		int
			numscan,
			w, x, y, z, c;

		std::string
			line;

		while (getline(in, line, '\n')) {
			line = tnstring::trim(line);
			if (line[0] != '#' && line != "") {
				if (firstline) {
					numscan = sscanf(line.c_str(), "%d\t%d\t%d\t%d",&w,&x,&y,&z);
					mtx = new matrix(4, w, x, y, z);
					firstline = false;
				} else {
					numscan = sscanf(line.c_str(), "%d\t%d\t%d\t%d\t%d",&w,&x,&y,&z,&c);
					mtx->put(w,x,y,z,(double)c);
				}
			}
		}
		return *mtx;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	matrix upload_rna_potential_matrix(std::istream &in) {
		matrix
			*mtx;

		bool
			firstline = true;

		int
			numscan,
			w, x, y, z;
		double
			c;

		std::string
			line;

		while (getline(in, line, '\n')) {
			line = tnstring::trim(line);
			if (line[0] != '#' && line != "") {
				if (firstline) {
					numscan = sscanf(line.c_str(), "%d\t%d\t%d\t%d",&w,&x,&y,&z);
					mtx = new matrix(4, w, x, y, z);
					firstline = false;
				} else {
					numscan = sscanf(line.c_str(), "%d\t%d\t%d\t%d\t%lf",&w,&x,&y,&z,&c);
					mtx->put(w,x,y,z,c);
				}
			}
		}
		return *mtx;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	void get_nrg_n_contacts(double *nrg, int *n_contacts, tnpdb::pdb in, matrix *pot) {
		matrix
			*mtx;
			
		std::vector<tnpot::point_type>
			nu_point_type;
	
		int
			n_types = 0,
			n_dist_classes = 0,
			n;
		
		n_types = tnpot::number_of_atom_types_per_label(tnpdb::residue::NU);
		n_dist_classes = tnpot::number_of_distance_classes();
		nu_point_type = get_point_types(in, tnpdb::residue::NU);
		
		std::vector<int>
			dims = pot->get_dimensions();
			
		if (dims[1] != n_types || dims[2] != n_types) {
			std::cerr << "ERROR: The number of atom types does not match the potential matrix." << std::endl;
			exit(1);	
		}	
			
		if (dims[3] != n_dist_classes) {
			std::cerr << "ERROR: The number of distance classes does not match the potential matrix." << std::endl;
			exit(1);
		}
		
		
		//std::cerr << dims[0] <<  " " << dims[1] <<  " " << dims[2] <<  " " << dims[3] <<  " -> " << n_types << " "<< n_dist_classes << std::endl; 
		mtx = new matrix(4, dims[0], dims[1], dims[2], dims[3]);
		get_rna_contact_matrix(mtx, &nu_point_type);
		
		*nrg = 0.0;
		*n_contacts = 0;

		for (unsigned int k = 0; k < dims[0]; ++k) {
			for (unsigned int i = 0; i < dims[1]; ++i) {
				for (unsigned int j = 0; j < dims[2]; ++j) {
					for (unsigned int d = 0; d < dims[3]; ++d) {
						n = (int)mtx->get(k, i, j, d);
						*nrg = *nrg +  n * pot->get(k, i, j, d);
						*n_contacts = *n_contacts + n;
					}
				}
			}
		}
		
		delete mtx, nu_point_type;
		
		return;
		
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	void randomize_resname_point_types(std::vector<point_type> *pt) {
		srand(time(0));
		
		unsigned int
			n_pt = pt->size(),
			rnd_num;
			
		std::string
			resname;
			
		std::vector<point_type>::iterator i, j;

		
		for (i = pt->begin(); i != pt->end(); ++i) {
			rnd_num = randomNum(0, n_pt - 1);
			resname = i->_resname;
			j = pt->begin() + rnd_num;
			i->_resname = j->_resname;
			j->_resname = resname;
		}
		
		return;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void reassign_types_point_types(std::vector<point_type> *pt) {
		
		unsigned int
			n_pt = pt->size(),
			rnd_num;
		int
			ind;
			
		std::string
			resname;
			
		std::vector<point_type>::iterator it;
		
		
		tnpdb::atom
			atom;
		
		for (it = pt->begin(); it != pt->end(); ++it) {
			
			atom.set_name(tnpdb::atom_string_to_name(it->_atmname, tnpdb::residue::NU));
			atom.set_ter(it->_ter);
			
			ind = get_atom_types_index(atom, tnpdb::residue_string_to_name(it->_resname));
			if (ind < 0) {
				std::cerr << "WARNING: In reassignment of types, cannot find new type definition. Setting the type = 1." << std::endl;
				it->_type = 1;
			} else {
				it->_type = get_atom_types_from_index(ind)[0];
			}
			
		}
		
		return;
	}

	
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	void get_nrg_mean_nrg_sd(double *nrg_mean, double *nrg_sd, tnpdb::pdb in, matrix *pot, int nz) {
		*nrg_mean = 0.0;
		*nrg_sd = 0.0;

		if (nz == 0) return;
		
		
		matrix
			*mtx;
			
		std::vector<tnpot::point_type>
			nu_point_type;
	
		int
			n_types = 0,
			n_dist_classes = 0,
			n;
		
		double
			cur_nrg;
			
		std::vector<double>
			nrg;
		
		n_types = tnpot::number_of_atom_types_per_label(tnpdb::residue::NU);
		n_dist_classes = tnpot::number_of_distance_classes();
		nu_point_type = get_point_types(in, tnpdb::residue::NU);
		
		std::vector<int>
			dims = pot->get_dimensions();
			
		if (dims[1] != n_types || dims[2] != n_types) {
			std::cerr << "ERROR: The number of atom types does not match the potential matrix." << std::endl;
			exit(1);	
		}	
			
		if (dims[3] != n_dist_classes) {
			std::cerr << "ERROR: The number of distance classes does not match the potential matrix." << std::endl;
			exit(1);
		}
		
		for (unsigned int z = 0; z < nz; z++) { 
			mtx = new matrix(4, dims[0], dims[1], dims[2], dims[3]);
			randomize_resname_point_types(&nu_point_type);
			reassign_types_point_types(&nu_point_type);
			
			
			//for(unsigned int p=0;p<nu_point_type.size();p++) {
			//	std::cerr << tnstring::trim(nu_point_type[p]._resname); 
			//}
			//std::cerr << std::endl;
			
			//for(unsigned int p=0;p<nu_point_type.size();p++) {
			//	std::cerr << nu_point_type[p]._type; 
			//}
			//std::cerr << std::endl;
			
			
			get_rna_contact_matrix(mtx, &nu_point_type);
			
			cur_nrg = 0.0;
			for (unsigned int k = 0; k < dims[0]; ++k) {
				for (unsigned int i = 0; i < dims[1]; ++i) {
					for (unsigned int j = 0; j < dims[2]; ++j) {
						for (unsigned int d = 0; d < dims[3]; ++d) {
							n = (int)mtx->get(k, i, j, d);
							cur_nrg += n * pot->get(k, i, j, d);
						}
					}
				}
			}
			
			nrg.push_back(cur_nrg);
			*nrg_mean = *nrg_mean + cur_nrg;
			
			//std::cerr << cur_nrg << std::endl;
		
			delete mtx;
		}
		
		*nrg_mean = *nrg_mean / nz;
		
		for (unsigned int i = 0; i < nz; i++) {
			*nrg_sd = *nrg_sd + (nrg[i] - *nrg_mean)*(nrg[i] - *nrg_mean);		
		}
		*nrg_sd = *nrg_sd / nz;
		*nrg_sd = sqrt(*nrg_sd);
		
		return;
	}
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<pick> get_nrg_profile(tnpdb::pdb in, matrix *pot, int w) {
						
		std::vector<tnpot::point_type>
			pt = get_point_types(in, tnpdb::residue::NU);
	
		float
			d = -1.0,
			val;
			
		int
			i, j,
			size = pt.size(),
			d_class = -1,
			k = 0,
			k_max = pot->get_dimensions()[0] - 1,
			min_resindex = 999999,
			max_resindex = -999999,
			first_resindex = 0,
			total = 0;
			
		std::vector<int>
			dims = pot->get_dimensions();
		
		for (i = 0; i < size; i++) {
			if (pt[i]._resindex < min_resindex) min_resindex = pt[i]._resindex;
			if (pt[i]._resindex > max_resindex) max_resindex = pt[i]._resindex;
			for (j = i + 1; j < size; j++) {
				
				d = tnpdb::distance(pt[i]._coords, pt[j]._coords);
				d_class = get_distance_class(d);
				
				if (d_class > 0) {
					
					k = pt[j]._resindex - pt[i]._resindex - 1;	
					
					if (k >= 0 && k < k_max) {
					
						val = pot->get(k, pt[i]._type - 1, pt[j]._type - 1, d_class - 1);
						pt[i]._value += val;
						pt[i]._value2 += 1.0;
						pt[j]._value += val;
						pt[j]._value2 += 1.0;
						total++;
					
					} else if (k >= k_max){
						
						val = pot->get(k_max, pt[i]._type - 1, pt[j]._type - 1, d_class - 1);
						pt[i]._value += val;
						pt[i]._value2 += 1.0;
						pt[j]._value += val;
						pt[j]._value2 += 1.0;
						total++;
					
					}
				}
			}
		}
		
		/*
		std::vector<float>
			profile(max_resindex - min_resindex + 1, 0.0),
			w_profile(max_resindex - min_resindex + 1, 0.0);
		
			
		std::vector<int>
			n_int(profile.size(),0);
		*/
			
		int
			cur_resindex = -9999;
		
		std::vector<pick>
			profile,
			w_profile;
			
		pick
			*p;

		//std::cerr << size << std::endl; 
			
		for (i = 0; i < size; i++) {
			if (pt[i]._resindex != cur_resindex) {
				p = new pick;
				p->_resindex = pt[i]._resindex;
				p->_resname = pt[i]._resname;
				p->_chn = pt[i]._original_chn;
				
				w_profile.push_back(*p);
				
				//p->_value = pt[i]._value/pt[i]._value2;
				p->_value = pt[i]._value;
				
				profile.push_back(*p);
				
				cur_resindex = pt[i]._resindex;
				
			} else {
				
				//profile.back()._value += pt[i]._value/pt[i]._value2;
				profile.back()._value += pt[i]._value;
			
			}
			
			//profile[pt[i]._resindex - min_resindex] += pt[i]._value/pt[i]._value2;
			//n_int[pt[i]._resindex - min_resindex] += 1;
			
		}
		
		/*
		for (i = 0; i < profile.size(); i++) {
			//std::cerr << "----> "<<profile[i] << "\t" << n_int[i]<<std::endl;
			profile[i] /= total;
			//std::cerr << "----> "<<profile[i] << std::endl;
		}
		*/
		
		/*
		size = profile.size();
		for (i = (w-1)/2; i < size - (w-1)/2; i++) {
			for (j = i - (w-1)/2; j <= i + (w-1)/2; j++) {
				w_profile[i] += profile[j];
			}
			w_profile[i] /= w;
		}
		for (i = 0; i < (w-1)/2; i++) {
			w_profile[i] = w_profile[(w-1)/2];
		}
		for (i = size - (w-1)/2; i < size; i++) {
			w_profile[i] = w_profile[size - (w-1)/2 - 1];
		}
		*/		

		size = profile.size();
		for (i = (w-1)/2; i < size - (w-1)/2; i++) {
			for (j = i - (w-1)/2; j <= i + (w-1)/2; j++) {
				w_profile[i]._value += profile[j]._value;
			}
			w_profile[i]._value /= w;
		}
		for (i = 0; i < (w-1)/2; i++) {
			w_profile[i]._value = w_profile[(w-1)/2]._value;
		}
		for (i = size - (w-1)/2; i < size; i++) {
			w_profile[i]._value = w_profile[size - (w-1)/2 - 1]._value;
		}		
		
		
		return w_profile;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	void asgl_potential_plot(matrix *pot) {

		std::vector<int>
			dims = pot->get_dimensions();
			
		std::string
			dir = "./plot_data_asgl/",
			subdir = "_data/",
			path;
		
		char
			filename[200],
			top_filename[50];
			
		int
			status = 0,
			numprint = 0,
			n_plot = 0,
			n_page = 1,
			n_top = 1;
			
		float
			min = FLAG_UP,
			max = FLAG_DOWN,
			dval;
			
		bool
			new_top = true;
			
		std::ofstream
			out,
			out_top;
			
			
					
		status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		path = dir + subdir;
		status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		
		
		for (unsigned int k = 0; k < dims[0]; ++k) {
			for (unsigned int i = 0; i < dims[1]; ++i) {
				for (unsigned int j = 0; j < dims[2]; ++j) {
					numprint = sprintf(filename, "%s%s_k%d_i%d_j%d.dat", dir.c_str(), subdir.c_str(), k, i, j);
					out.open(filename, std::iostream::out);
					
					for (unsigned int d = 0; d < dims[3]; ++d) {
						dval = pot->get(k,i,j,d);
						out << (d + 0.5) << "\t" << std::setprecision(4) << std::setiosflags(std::ios::fixed) << dval << std::endl;
						if (min > dval) min = dval;
						if (max < dval) max = dval;
					}
					
					out.close();
					
				}
			}
		}

		numprint = sprintf(filename, "%s%s_zero.dat", dir.c_str(), subdir.c_str());
		out.open(filename, std::iostream::out);
		out << "0\t0.0" << std::endl;
		out << dims[3] << "\t0.0" << std::endl;
		out.close();


		for (unsigned int k = 0; k < dims[0]; ++k) {
			for (unsigned int i = 0; i < dims[1]; ++i) {
				for (unsigned int j = 0; j < dims[2]; ++j) {
					
					if (new_top) {
						numprint = sprintf(top_filename, "_%d.top", n_top);
						numprint = sprintf(filename, "%s%s", dir.c_str(), top_filename);
						out_top.open(filename, std::iostream::out);
						if (out_top.fail() || !out_top.is_open()) {
							std::cerr << "ERROR: Cannot open file \"" << filename << "\" to write." << std::endl;
							exit(1);
						}
						new_top = false;
					}
					
					out_top << "RESET\n";
					
					numprint = sprintf(filename, "./%s_k%d_i%d_j%d.dat", subdir.c_str(), k, i, j);
					
					out_top << "READ_TABLE FILE '" << filename << "'" << std::endl;
					out_top << "SET POSITION " << 16 + n_plot << " 0" << std::endl;
					out_top << "SET TICK_FONT 4" << std::endl;
					out_top << "SET X_TICK_DECIMALS -1" << std::endl;
					out_top << "SET WORLD_WINDOW 0 "; 
					out_top << std::setprecision(2) << std::setiosflags(std::ios::fixed) << min - 0.01 << " ";
					out_top << std::setprecision(0) << std::setiosflags(std::ios::fixed) << dims[3] << " ";
					out_top << std::setprecision(2) << std::setiosflags(std::ios::fixed) << max + 0.01 << std::endl;
					out_top << "WORLD" << std::endl;
					out_top << "AXES2D" << std::endl;
					out_top << "RESET_CAPTIONS" << std::endl;
					out_top << "CAPTION CAPTION_POSITION 1, CAPTION_FONT 4, CAPTION_TEXT " << "'k = " << k << "  i = " << i << "  j = " << j << "'" << std::endl;
					out_top << "CAPTION CAPTION_POSITION 2, CAPTION_FONT 4, CAPTION_TEXT 'distance (Ang)'" << std::endl;
					out_top << "CAPTION CAPTION_POSITION 3, CAPTION_FONT 4, CAPTION_TEXT '@D@E (Kcal/mol)'" << std::endl;
					out_top << "RESET_LEGEND" << std::endl;
					out_top << "SET PLOT2D_LINE_TYPE = 1" << std::endl;
					out_top << "PLOT2D" << std::endl;
					
					numprint = sprintf(filename, "./%s_zero.dat", subdir.c_str());
					
					out_top << "READ_TABLE FILE '" << filename << "'" << std::endl;
					out_top << "SET PLOT2D_LINE_TYPE = 2" << std::endl;
					out_top << "PLOT2D" << std::endl;
					
					n_plot++;
					
					if (n_plot == 32) {
					
						n_plot = 0;
						n_page++;
						
						if (n_page == 15) {
						
							n_page = 1;
							new_top = true;
							n_top++;
							out_top.close();
						
						} else {
						
							out_top << "NEW_PAGE" << std::endl;
						}	
					}
					
					if (new_top) {
						
						numprint = sprintf(filename, "%s%s", dir.c_str(), "execute.sh");
						out.open(filename, std::iostream::app | std::iostream::ate);
						if (out.fail() || !out.is_open()) {
							std::cerr << "ERROR: Cannot open file \"" << filename << "\" to write." << std::endl;
							exit(1);
						}
						status = chmod(filename, S_IRWXU | S_IRWXG | S_IRWXO);
						
						out << "asgl " << top_filename << std::endl;
						
						out.close();
					}
				}
			}
		}
		

		return;
	}	

	//////////////////////////////////////////////////////////////////////////////////////////////
	void asgl_potential_plots_by_k(std::vector<matrix> *pots, std::vector<int> ks) {

		std::vector<int>
			dims;
			
		std::string
			dir = "./plot_data_asgl_by_k/",
			subdir = "_data/",
			path;
		
		char
			filename[200],
			top_filename[50];
			
		int
			status = 0,
			numprint = 0,
			n_plot = 0,
			n_page = 1,
			n_top = 1,
			k;
			
		float
			min = FLAG_UP,
			max = FLAG_DOWN,
			dval;
			
		bool
			new_top = true;
			
		std::ofstream
			out,
			out_top;
		
		
		std::vector<matrix>::iterator it;	
			
					
		status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		path = dir + subdir;
		status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		
		k = 0;
		
		for (it = pots->begin(); it != pots->end(); it++) {
			dims = it->get_dimensions();
			
			if (dims[0] == ks[k] + 1) {
				for (unsigned int i = 0; i < dims[1]; ++i) {
					for (unsigned int j = 0; j < dims[2]; ++j) {
						numprint = sprintf(filename, "%s%s_k%d_i%d_j%d.dat", dir.c_str(), subdir.c_str(), ks[k], i, j);
						out.open(filename, std::iostream::out);
				
						for (unsigned int d = 0; d < dims[3]; ++d) {
							dval = it->get(ks[k],i,j,d);
							out << d << "\t" << std::setprecision(4) << std::setiosflags(std::ios::fixed) << dval << std::endl;
							if (min > dval) min = dval;
							if (max < dval) max = dval;
						}
				
						out.close();	
					}
				}
			} else {
				std::cerr << "ERROR: " << std::endl;
				exit(1);
			}
			k++;
		}

		numprint = sprintf(filename, "%s%s_zero.dat", dir.c_str(), subdir.c_str());
		out.open(filename, std::iostream::out);
		out << "0\t0.0" << std::endl;
		out << dims[3] << "\t0.0" << std::endl;
		out.close();


		k = 0;
		for (it = pots->begin(); it != pots->end(); it++) {
			dims = it->get_dimensions();
			if (dims[0] == ks[k] + 1) {
				for (unsigned int i = 0; i < dims[1]; ++i) {
					for (unsigned int j = 0; j < dims[2]; ++j) {
					
						if (new_top) {
							numprint = sprintf(top_filename, "_%d.top", n_top);
							numprint = sprintf(filename, "%s%s", dir.c_str(), top_filename);
							out_top.open(filename, std::iostream::out);
							if (out_top.fail() || !out_top.is_open()) {
								std::cerr << "ERROR: Cannot open file \"" << filename << "\" to write." << std::endl;
								exit(1);
							}
							new_top = false;
						}
					
						out_top << "RESET\n";
					
						numprint = sprintf(filename, "./%s_k%d_i%d_j%d.dat", subdir.c_str(), ks[k], i, j);
					
						out_top << "READ_TABLE FILE '" << filename << "'" << std::endl;
						out_top << "SET POSITION " << 16 + n_plot << " 0" << std::endl;
						out_top << "SET TICK_FONT 4" << std::endl;
						out_top << "SET X_TICK_DECIMALS -1" << std::endl;
						out_top << "SET WORLD_WINDOW 0 "; 
						out_top << std::setprecision(2) << std::setiosflags(std::ios::fixed) << min - 0.01 << " ";
						out_top << std::setprecision(0) << std::setiosflags(std::ios::fixed) << dims[3] << " ";
						out_top << std::setprecision(2) << std::setiosflags(std::ios::fixed) << max + 0.01 << std::endl;
						out_top << "WORLD" << std::endl;
						out_top << "AXES2D" << std::endl;
						out_top << "RESET_CAPTIONS" << std::endl;
						out_top << "CAPTION CAPTION_POSITION 1, CAPTION_FONT 4, CAPTION_TEXT " << "'k >= " << ks[k] << "  i = " << i << "  j = " << j << "'" << std::endl;
						out_top << "CAPTION CAPTION_POSITION 2, CAPTION_FONT 4, CAPTION_TEXT 'distance (Ang)'" << std::endl;
						out_top << "CAPTION CAPTION_POSITION 3, CAPTION_FONT 4, CAPTION_TEXT '@D@E (Kcal/mol)'" << std::endl;
						out_top << "RESET_LEGEND" << std::endl;
						out_top << "SET PLOT2D_LINE_TYPE = 18" << std::endl;
						out_top << "PLOT2D" << std::endl;
					
						numprint = sprintf(filename, "./%s_zero.dat", subdir.c_str());
					
						out_top << "READ_TABLE FILE '" << filename << "'" << std::endl;
						out_top << "SET PLOT2D_LINE_TYPE = 19" << std::endl;
						out_top << "PLOT2D" << std::endl;
					
						n_plot++;
					
						if (n_plot == 32) {
					
							n_plot = 0;
							n_page++;
						
							if (n_page == 15) {
						
								n_page = 1;
								new_top = true;
								n_top++;
								out_top.close();
						
							} else {
							
								out_top << "NEW_PAGE" << std::endl;
							}	
						}
					
						if (new_top) {
						
							numprint = sprintf(filename, "%s%s", dir.c_str(), "execute.sh");
							out.open(filename, std::iostream::app | std::iostream::ate);
							if (out.fail() || !out.is_open()) {
								std::cerr << "ERROR: Cannot open file \"" << filename << "\" to write." << std::endl;
								exit(1);
							}
							status = chmod(filename, S_IRWXU | S_IRWXG | S_IRWXO);
							
							out << "asgl " << top_filename << std::endl;
							
							out.close();
						}
					}
				}
			} else {
				std::cerr << "ERROR: " << std::endl;
				exit(1);
			}
			k++;
		}
		
		return;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	void print_rasmol_script(std::vector<pick> *profile) {
		
		std::vector<pick>::iterator
			it;
		
		
		std::ofstream
			out;
			
		std::vector<rgb>
			stops;
		
		
		float
			min_nrg = FLAG_UP,
			max_nrg = FLAG_DOWN;
			
		unsigned int
			size = profile->size();
			
		if (size == 0) {
			std::cerr << "ERROR: Not profile set yet." << std::endl;
			exit(1);
		}
		
		for (it = profile->begin(); it != profile->end(); it++) {
			if (it->_value < min_nrg) min_nrg = it->_value;
			if (it->_value > max_nrg) max_nrg = it->_value;
		}
		
		stops.push_back((rgb){0,0,255}); // blue
		stops.push_back((rgb){255,255,255}); // white
		stops.push_back((rgb){255,0,0}); // red
		
		Gradient grad(min_nrg, max_nrg, stops);
		
		
		rgb
			color;
			
		out.open("profile.scr", std::iostream::out);
		
		
		out << "echo Min energy (blue): " << min_nrg << std::endl;
		out << "echo Max energy (red) : " << max_nrg << std::endl;
		out << "wireframe 90" << std::endl;
		
		
		for (it = profile->begin(); it != profile->end(); it++) {
			out << "select " << it->_resindex << std::endl;
			color = grad.getRgb(it->_value);
			out << "color [" << color.red << "," << color.green << "," << color.blue << "]" << std::endl;
		}
		out.close();
		
	}
	
};



