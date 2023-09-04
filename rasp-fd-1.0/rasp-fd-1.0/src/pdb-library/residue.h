#include <vector>
#include <map>
#include "atom.h"

namespace tnpdb {

	class residue {
		friend class pdb;
		friend class chain;
		public:
			enum label {T_INV, AA, NU, WA, ME};
			enum name {
				N_INV, AAB, NUB, //2
				// Aminoacids, Centroid
				ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, CEA, //23
				// Nucleotides, ATP, Centroid
				ADE, CYT, GUA, THY, URA, ANY, PUR, PYR, ATP, CEN, //33
				// Metals
				MG, ZN, CA, NA, NI, CL, FE, //40
				// Water
				HOH //41
			};

			typedef std::pair< index<atom>, index<atom> > bond;

			//Constructor
			residue(){}
			residue(name);
			//Retrieve data
			name get_name() const;                        // The name for the residue
			label get_label() const;
			char get_chain() const;
			unsigned int get_number_of_atoms() const;     // The number of atoms present in the residue
			const atom &get_atom(atom::name) const;       // Return the data for an atom
			const atom &get_atom_at(int) const;
			atom &get_atom_at(int);
			atom &get_atom_to_modify(atom::name);
			index<atom> get_last_atom_index() const;      // The index of the last atom in the residue
			index<residue> get_index() const;             // The index for the residue
			//Set data
			void set_index(index<residue>);               // Set the index for the residue
			void set_label(label);
			void set_chain(char);
			void add_atom(const atom &);
			//Output
			void write_pdb(char, std::ostream &) const;   // Write the lines for a pdb file

		protected: // Accessible by members and childs
			index<atom> get_min_atom_index() const;
			bool can_have_atom(atom::name) const;         // Return true if the residue can have atoms of this name
		private: // Accessible by members and friends
			std::vector<atom>
				_atoms;
			
			std::vector<bond>
				_bonds;
				
			name
				_name;

			label
				_label;
			
			index<residue>
				_index;
				
			index<atom>
				_min_atom_index;

			char
				_chain;
	};

	// Conversion data
	atom::etype atom_name_to_etype(atom::name, residue::label);
	std::string atom_name_to_string(atom::name, residue::label);
	atom::name atom_string_to_name(const std::string &, residue::label);
	residue::name residue_string_to_name(const std::string &);
	std::string residue_name_to_string(residue::name);
	residue::label residue_name_to_label(residue::name);
	std::string residue_name_to_1_letter_string(residue::name);

	//Other
	atom_data *get_atom_name_data(residue::label);
	void initialize_hvatm_res_data();

    //Comparison class
	class atomp_index_less_than{
		public:
			inline bool operator()(const std::pair< atom::name, atom > &a1, const std::pair< atom::name, atom > &a2){
				return (a1.second.get_index() < a2.second.get_index());
			}
	};


	// RESIDUE DATA
    struct _res_data {
		char s[4];
		residue::name name;
		residue::label label;
	};

	typedef _res_data res_data;

	static res_data res_name_data[] = {
		// Aminoacids
		{"ALA", residue::ALA, residue::AA},
		{"ARG", residue::ARG, residue::AA},
		{"ASN", residue::ASN, residue::AA},
		{"ASP", residue::ASP, residue::AA},
		{"CYS", residue::CYS, residue::AA},
		{"GLN", residue::GLN, residue::AA},
		{"GLU", residue::GLU, residue::AA},
		{"GLY", residue::GLY, residue::AA},
		{"HIS", residue::HIS, residue::AA},
		{"ILE", residue::ILE, residue::AA},
		{"LEU", residue::LEU, residue::AA},
		{"LYS", residue::LYS, residue::AA},
		{"MET", residue::MET, residue::AA},
		{"PHE", residue::PHE, residue::AA},
		{"PRO", residue::PRO, residue::AA},
		{"SER", residue::SER, residue::AA},
		{"THR", residue::THR, residue::AA},
		{"TRP", residue::TRP, residue::AA},
		{"TYR", residue::TYR, residue::AA},
		{"VAL", residue::VAL, residue::AA},
		{"CEA", residue::CEA, residue::AA},
		// Nucleotides
		{"  A", residue::ADE, residue::NU}, {" DA", residue::ADE, residue::NU},
		{"  C", residue::CYT, residue::NU}, {" DC", residue::CYT, residue::NU},
		{"  G", residue::GUA, residue::NU}, {" DG", residue::GUA, residue::NU},
		{"  T", residue::THY, residue::NU}, {" DT", residue::THY, residue::NU},
		{"  U", residue::URA, residue::NU},
		{"  N", residue::ANY, residue::NU},
		{"  R", residue::PUR, residue::NU},
		{"  Y", residue::PYR, residue::NU},
		{"CEN", residue::CEN, residue::NU},
		// Metals
		{" MG", residue::MG, residue::ME},
		{" ZN", residue::ZN, residue::ME},
		{" CA", residue::CA, residue::ME},
		{" NA", residue::NA, residue::ME},
		{" NI", residue::NI, residue::ME},
		{" CL", residue::CL, residue::ME},
		{" FE", residue::FE, residue::ME},
		// Water
		{"HOH", residue::HOH, residue::WA},
		// ATP
		//{"ATP", residue::ATP, residue::NU},
		// Other
		{"UKN", residue::N_INV, residue::T_INV}

	};

	static res_data res_1_letter_name_data[] = {
		// Aminoacids
		{"A", residue::ALA, residue::AA},
		{"R", residue::ARG, residue::AA},
		{"N", residue::ASN, residue::AA},
		{"D", residue::ASP, residue::AA},
		{"C", residue::CYS, residue::AA},
		{"Q", residue::GLN, residue::AA},
		{"E", residue::GLU, residue::AA},
		{"G", residue::GLY, residue::AA},
		{"H", residue::HIS, residue::AA},
		{"I", residue::ILE, residue::AA},
		{"L", residue::LEU, residue::AA},
		{"K", residue::LYS, residue::AA},
		{"M", residue::MET, residue::AA},
		{"F", residue::PHE, residue::AA},
		{"P", residue::PRO, residue::AA},
		{"S", residue::SER, residue::AA},
		{"T", residue::THR, residue::AA},
		{"W", residue::TRP, residue::AA},
		{"Y", residue::TYR, residue::AA},
		{"V", residue::VAL, residue::AA},
		// Nucleotides
		{"A", residue::ADE, residue::NU},
		{"C", residue::CYT, residue::NU},
		{"G", residue::GUA, residue::NU},
		{"T", residue::THY, residue::NU},
		{"U", residue::URA, residue::NU},
		{"N", residue::ANY, residue::NU},
		{"R", residue::PUR, residue::NU},
		{"Y", residue::PYR, residue::NU},
		{"-", residue::N_INV, residue::T_INV}
	};


	static std::map<residue::name, atom::name *> hvatm_res_data;
	//static std::map<residue::name, atom::name *> hyatm_res_data;

	// DUMMY
	static atom::name hvatm_dum[] =  {atom::N_INV};

	//AA BACKBONE
	static atom::name hvatm_aab[] =  {atom::AA_N, atom::AA_CA, atom::AA_C, atom::AA_O, atom::AA_OXT, atom::N_INV};
	static atom::name hyatm_aab[] =  {atom::AA_H1, atom::AA_H2, atom::AA_H3, atom::AA_H, atom::AA_HA, atom::AA_HA2, atom::AA_HA3, atom::N_INV};

	//ALA
	static atom::name hvatm_ala[] =  {atom::AA_CB, atom::N_INV};
	static atom::name hyatm_ala[] =  {atom::AA_HB1, atom::AA_HB2, atom::AA_HB3, atom::N_INV};

	//ARG
	static atom::name hvatm_arg[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD, atom::AA_NE, atom::AA_CZ, atom::AA_NH1, atom::AA_NH2, atom::N_INV};
	static atom::name hyatm_arg[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG2, atom::AA_HG3, atom::AA_HD2, atom::AA_HD3, atom::AA_HE, atom::AA_HH11, atom::AA_HH12, atom::AA_HH21, atom::AA_HH22, atom::N_INV};

	//ASN
	static atom::name hvatm_asn[] =  {atom::AA_CB, atom::AA_CG, atom::AA_OD1, atom::AA_ND2, atom::N_INV};
	static atom::name hyatm_asn[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HD21, atom::AA_HD22, atom::N_INV};

	//ASP
	static atom::name hvatm_asp[] =  {atom::AA_CB, atom::AA_CG, atom::AA_OD1, atom::AA_OD2, atom::N_INV};
	static atom::name hyatm_asp[] =  {atom::AA_HB2, atom::AA_HB3, atom::N_INV};

	//CYS
	static atom::name hvatm_cys[] =  {atom::AA_CB, atom::AA_SG, atom::N_INV};
	static atom::name hyatm_cys[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG, atom::N_INV};

	//GLN
	static atom::name hvatm_gln[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD, atom::AA_OE1, atom::AA_NE2, atom::N_INV};
	static atom::name hyatm_gln[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG2, atom::AA_HG3, atom::AA_HE21, atom::AA_HE22, atom::N_INV};

	//GLN
	static atom::name hvatm_glu[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD, atom::AA_OE1, atom::AA_OE2, atom::N_INV};
	static atom::name hyatm_glu[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG2, atom::AA_HG3, atom::N_INV};

	//GLY
	static atom::name hvatm_gly[] =  {atom::N_INV};
	static atom::name hyatm_gly[] =  {atom::AA_HA2, atom::AA_HA3, atom::N_INV};

	//HIS
	static atom::name hvatm_his[] =  {atom::AA_CB, atom::AA_CG, atom::AA_ND1, atom::AA_CD2, atom::AA_CE1, atom::AA_NE2, atom::N_INV};
	static atom::name hyatm_his[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HD1, atom::AA_HD2, atom::AA_HE1, atom::N_INV};

	//ILE
	static atom::name hvatm_ile[] =  {atom::AA_CB, atom::AA_CG1, atom::AA_CG2, atom::AA_CD1, atom::N_INV};
	static atom::name hyatm_ile[] =  {atom::AA_HG12, atom::AA_HG13, atom::AA_HG21, atom::AA_HG22, atom::AA_HG23, atom::AA_HD11, atom::AA_HD12, atom::AA_HD13, atom::N_INV};

	//LEU
	static atom::name hvatm_leu[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD1, atom::AA_CD2, atom::N_INV};
	static atom::name hyatm_leu[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG, atom::AA_HD11, atom::AA_HD12, atom::AA_HD13, atom::AA_HD21, atom::AA_HD22, atom::AA_HD23, atom::N_INV};

	//LYS
	static atom::name hvatm_lys[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD, atom::AA_CE, atom::AA_NZ, atom::N_INV};
	static atom::name hyatm_lys[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG2, atom::AA_HG3, atom::AA_HD2, atom::AA_HD3, atom::AA_HE2, atom::AA_HE3, atom::AA_HZ1, atom::AA_HZ2, atom::AA_HZ3, atom::N_INV};

	//MET
	static atom::name hvatm_met[] =  {atom::AA_CB, atom::AA_CG, atom::AA_SD, atom::AA_CE, atom::N_INV};
	static atom::name hyatm_met[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG2, atom::AA_HG3, atom::AA_HE1, atom::AA_HE2, atom::AA_HE3, atom::N_INV};

	//PHE
	static atom::name hvatm_phe[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD1, atom::AA_CD2, atom::AA_CE1, atom::AA_CE2, atom::AA_CZ, atom::N_INV};
	static atom::name hyatm_phe[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HD1, atom::AA_HD2, atom::AA_HE1, atom::AA_HE2, atom::AA_HZ, atom::N_INV};

	//PRO
	static atom::name hvatm_pro[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD, atom::N_INV};
	static atom::name hyatm_pro[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG2, atom::AA_HG3, atom::AA_HD2, atom::AA_HD3, atom::N_INV};

	//SER
	static atom::name hvatm_ser[] =  {atom::AA_CB, atom::AA_OG, atom::N_INV};
	static atom::name hyatm_ser[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HG, atom::N_INV};

	//THR
	static atom::name hvatm_thr[] =  {atom::AA_CB, atom::AA_OG1, atom::AA_CG2, atom::N_INV};
	static atom::name hyatm_thr[] =  {atom::AA_HB, atom::AA_HG1, atom::AA_HG21, atom::AA_HG22, atom::AA_HG23, atom::N_INV};

	//TRP
	static atom::name hvatm_trp[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD1, atom::AA_CD2, atom::AA_NE1, atom::AA_CE2, atom::AA_CE3, atom::AA_CZ2, atom::AA_CZ3, atom::AA_CH2, atom::N_INV};
	static atom::name hyatm_trp[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HD1, atom::AA_HE1, atom::AA_HE3, atom::AA_HZ2, atom::AA_HZ3, atom::AA_HH2, atom::N_INV};

	//TYR
	static atom::name hvatm_tyr[] =  {atom::AA_CB, atom::AA_CG, atom::AA_CD1, atom::AA_CD2, atom::AA_CE1, atom::AA_CE2, atom::AA_CZ, atom::AA_OH, atom::N_INV};
	static atom::name hyatm_tyr[] =  {atom::AA_HB2, atom::AA_HB3, atom::AA_HD1, atom::AA_HD2, atom::AA_HE1, atom::AA_HE2, atom::AA_HH, atom::N_INV};

	//VAL
	static atom::name hvatm_val[] =  {atom::AA_CB, atom::AA_CG1, atom::AA_CG2, atom::N_INV};
	static atom::name hyatm_val[] =  {atom::AA_HB, atom::AA_HG11, atom::AA_HG12, atom::AA_HG13, atom::AA_HG21, atom::AA_HG22, atom::AA_HG23, atom::N_INV};

	//NU BACKBONE
	static atom::name hvatm_nub[] =  {atom::NU_P, atom::NU_OP1, atom::NU_OP2, atom::NU_OP3, atom::NU_O5p, atom::NU_C5p, atom::NU_C4p, atom::NU_O4p, atom::NU_C3p, atom::NU_O3p, atom::NU_C2p, atom::NU_O2p, atom::NU_C1p, atom::N_INV};
	static atom::name hyatm_nub[] =  {atom::NU_H5p, atom::NU_H5pp, atom::NU_H4p, atom::NU_H3p, atom::NU_HO3p, atom::NU_H2p, atom::NU_H2pp, atom::NU_HO2p, atom::NU_H1p, atom::NU_HO5p, atom::N_INV};
	
	//ADE
	static atom::name hvatm_ade[] =  {atom::NU_N9, atom::NU_C8, atom::NU_N7, atom::NU_C5, atom::NU_C6, atom::NU_N6, atom::NU_N1, atom::NU_C2, atom::NU_N3, atom::NU_C4, atom::N_INV};
	static atom::name hyatm_ade[] =  {atom::NU_H8, atom::NU_H61, atom::NU_H62, atom::NU_H2, atom::N_INV};

	//CYT
	static atom::name hvatm_cyt[] =  {atom::NU_N1, atom::NU_C2, atom::NU_O2, atom::NU_N3, atom::NU_C4, atom::NU_N4, atom::NU_C5, atom::NU_C6, atom::N_INV};
	static atom::name hyatm_cyt[] =  {atom::NU_H41, atom::NU_H42, atom::NU_H5, atom::NU_H6, atom::N_INV};

	//GUA
	static atom::name hvatm_gua[] =  {atom::NU_N9, atom::NU_C8, atom::NU_N7, atom::NU_C5, atom::NU_C6, atom::NU_O6, atom::NU_N1, atom::NU_C2, atom::NU_N2, atom::NU_N3, atom::NU_C4, atom::N_INV};
	static atom::name hyatm_gua[] =  {atom::NU_H8, atom::NU_H1, atom::NU_H21, atom::NU_H22, atom::N_INV};

	//THY
	static atom::name hvatm_thy[] =  {atom::NU_N1, atom::NU_C2, atom::NU_O2, atom::NU_N3, atom::NU_C4, atom::NU_O4, atom::NU_C5, atom::NU_C7, atom::NU_C6, atom::N_INV};
	static atom::name hyatm_thy[] =  {atom::NU_H3, atom::NU_H71, atom::NU_H72, atom::NU_H73, atom::NU_H6, atom::N_INV};

	//URA
	static atom::name hvatm_ura[] =  {atom::NU_N1, atom::NU_C2, atom::NU_O2, atom::NU_N3, atom::NU_C4, atom::NU_O4, atom::NU_C5, atom::NU_C6, atom::N_INV};
	static atom::name hyatm_ura[] =  {atom::NU_H3, atom::NU_H5, atom::NU_H6, atom::N_INV};

	//ATP
	static atom::name hvatm_atp[] =  {atom::NU_PA, atom::NU_PB, atom::NU_PG, atom::NU_O1A, atom::NU_O2A, atom::NU_O3A, atom::NU_O1B, atom::NU_O2B, atom::NU_O3B, atom::NU_O1G, atom::NU_O2G, atom::NU_O3G, atom::NU_O5p, atom::NU_C5p, atom::NU_C4p, atom::NU_O4p, atom::NU_C3p, atom::NU_O3p, atom::NU_C2p, atom::NU_O2p, atom::NU_C1p, atom::NU_N9, atom::NU_C8, atom::NU_N7, atom::NU_C5, atom::NU_C6, atom::NU_N6, atom::NU_N1, atom::NU_C2, atom::NU_N3, atom::NU_C4, atom::N_INV};
	//static atom::name hvatm_atp[] =  {atom::NU_PA, atom::NU_PB, atom::NU_PG, atom::NU_O1A, atom::NU_O2A, atom::NU_O3A, atom::NU_O1B, atom::NU_O2B, atom::NU_O3B, atom::NU_O1G, atom::NU_O2G, atom::NU_O3G, atom::NU_O5a, atom::NU_C5a, atom::NU_C4a, atom::NU_O4a, atom::NU_C3a, atom::NU_O3a, atom::NU_C2a, atom::NU_O2a, atom::NU_C1a, atom::NU_N9, atom::NU_C8, atom::NU_N7, atom::NU_C5, atom::NU_C6, atom::NU_N6, atom::NU_N1, atom::NU_C2, atom::NU_N3, atom::NU_C4, atom::N_INV};
	static atom::name hyatm_atp[] =  {atom::N_INV};

	//HOH
	static atom::name hvatm_hoh[] =  {atom::WA_O, atom::N_INV};
	static atom::name hyatm_hoh[] =  {atom::N_INV};
	
	//CEA
	static atom::name hvatm_cea[] =  {atom::AA_CEA, atom::N_INV};
	static atom::name hyatm_cea[] =  {atom::N_INV};

	//CEN
	static atom::name hvatm_cen[] =  {atom::NU_CEN, atom::N_INV};
	static atom::name hyatm_cen[] =  {atom::N_INV};
	
	//METALS
	static atom::name hvatm_mg[] =  {atom::ME_MG, atom::N_INV};
	static atom::name hyatm_mg[] =  {atom::N_INV};
	
	static atom::name hvatm_zn[] =  {atom::ME_ZN, atom::N_INV};
	static atom::name hyatm_zn[] =  {atom::N_INV};
	
	static atom::name hvatm_ca[] =  {atom::ME_CA, atom::N_INV};
	static atom::name hyatm_ca[] =  {atom::N_INV};
	
	static atom::name hvatm_na[] =  {atom::ME_NA, atom::N_INV};
	static atom::name hyatm_na[] =  {atom::N_INV};
	
	static atom::name hvatm_ni[] =  {atom::ME_NI, atom::N_INV};
	static atom::name hyatm_ni[] =  {atom::N_INV};
	
	static atom::name hvatm_cl[] =  {atom::ME_CL, atom::N_INV};
	static atom::name hyatm_cl[] =  {atom::N_INV};
	
	static atom::name hvatm_fe[] =  {atom::ME_FE, atom::N_INV};
	static atom::name hyatm_fe[] =  {atom::N_INV};


};

