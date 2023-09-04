#include <string>
#include <cassert>
#include <cstdio>
#include "point.h"
#include "index.h"

namespace tnpdb {

	class atom {
		friend class pdb;
		friend class residue;
		public:
			// The etype (element) of an atom
			enum etype {T_INV, C, N, H, O, S, P, MG, ZN, CA, NA, NI, CL, FE};
			enum name {
				N_INV,
				// Aminoacids
				// Carbons
				AA_C, AA_CA, AA_CB, AA_CD, AA_CD1, AA_CD2, AA_CE, AA_CE1, AA_CE2, AA_CE3, AA_CG, AA_CG1, AA_CG2, AA_CH2, AA_CZ, AA_CZ2, AA_CZ3,
				// Hydrogens
				AA_H, AA_H1, AA_H2, AA_H3, AA_HA, AA_HA2, AA_HA3, AA_HB, AA_HB1, AA_HB2, AA_HB3, AA_HD1, AA_HD11, AA_HD12, AA_HD13, AA_HD2,
				AA_HD21, AA_HD22, AA_HD23, AA_HD3, AA_HE, AA_HE1, AA_HE2, AA_HE21, AA_HE22, AA_HE3, AA_HG, AA_HG1, AA_HG11, AA_HG12, AA_HG13,
				AA_HG2, AA_HG21, AA_HG22, AA_HG23, AA_HG3, AA_HH, AA_HH11, AA_HH12, AA_HH2, AA_HH21, AA_HH22, AA_HZ, AA_HZ1, AA_HZ2, AA_HZ3,
				// Nitrogens
				AA_N, AA_ND1, AA_ND2, AA_NE, AA_NE1, AA_NE2, AA_NH1, AA_NH2, AA_NZ,
				// Oxygens
				AA_O, AA_OD1, AA_OD2, AA_OE1, AA_OE2, AA_OG, AA_OG1, AA_OH, AA_OXT, AA_CEA, //Centroid
				// Sulphurs
				AA_SD, AA_SG,

				// Nucleotides (p: prime, a: asterisc)
				// Carbonis
				NU_CA,
				NU_C1p, NU_C2, NU_C2p, NU_C3p, NU_C4, NU_C4p, NU_C5, NU_C5p, NU_C6, NU_C7, NU_C8,
				// Carbons ATP
				//NU_C1a, NU_C2a, NU_C3a, NU_C4a, NU_C5a,
				// Hydrogens
				NU_H1, NU_H1p, NU_H2, NU_H21, NU_H22, NU_H2p, NU_H2pp, NU_H3, NU_H3p, NU_H41, NU_H42, NU_H4p, NU_H5,
				NU_H5p, NU_H5pp, NU_H6, NU_H61, NU_H62, NU_H71, NU_H72, NU_H73, NU_H8, NU_HO2p, NU_HO3p, NU_HO5p,
				// Nitrogens
				NU_N1, NU_N2, NU_N3, NU_N4, NU_N6, NU_N7, NU_N9,
				// Oxygens
				NU_O2, NU_O2p, NU_O3p, NU_O4, NU_O4p, NU_O5p, NU_O6, NU_OP1, NU_OP2, NU_OP3, NU_CEN, //Centroid
				// Oxygens ATP
				NU_O1A, NU_O2A, NU_O3A, NU_O1B, NU_O2B, NU_O3B, NU_O1G, NU_O2G, NU_O3G, /*NU_O2a, NU_O3a, NU_O4a, NU_O5a,*/
				// Phosphorus
				NU_P,
				// Phosphorus ATP
				NU_PA, NU_PB, NU_PG,			
				 
				// Water General
				WA_CA,
				/*NU_O,*/ WA_O,
				
				//Metals
				ME_MG, ME_ZN, ME_CA, ME_NA, ME_NI, ME_CL, ME_FE
				
			};
			// Constructor (invalid atom)
			atom();
			// Overloadiing operator functions
			bool operator < (const atom &) const;
			bool operator == (const atom &) const;
			bool operator != (const atom &) const;
			// Methods for retrieving data
			index<atom> get_index() const;
			const point &get_coords() const;
			float get_occupancy() const;
			float get_temp_factor() const;
			const char *get_segment_id() const;
			const char *get_element() const;
			const char *get_charge() const;
			etype get_etype() const;
			name get_name() const;
			const char get_rec() const;
			const char get_alt_loc() const;
			bool  get_ter() const;
			float get_vdW_radii() const;
			// Methods for setting data
			void set_index(index<atom>);
			void set_coords(const point &);
			void set_occupancy(float);
			void set_temp_factor(float);
			void set_segment_id(const char *);
			void set_element(const char *);
			void set_charge(const char *);
			void set_etype(etype);
			void set_name(name);
			void set_rec(char);
			void set_alt_loc(char);
			void set_ter(bool);
		protected:
			// Members
			index<atom>
				_index;

			etype
				_etype;

			name
				_name;

			point
				_coords;

			float
				_occupancy,
				_temp_factor;

			std::string
				_segment_id,
				_element,
				_charge;

			char
				_alt_loc,
				_rec;

			bool
				_ter;
	};


	//Comparison class
	class atom_index_less_than{
		public:
			inline bool operator()(const atom &a1, const atom &a2){
				return (a1.get_index() < a2.get_index());
			}
	};

	// ATOM DATA
	struct _atom_data {
		char s[5];
		atom::name name;
		atom::etype etype;
	};

	typedef _atom_data atom_data;

	static atom_data atom_name_aa_data[] = {
		// Aminoacids
		{" C  ", atom::AA_C, atom::C},
		{" CA ", atom::AA_CA, atom::C},
		{" CB ", atom::AA_CB, atom::C},
		{" CD ", atom::AA_CD, atom::C},
		{" CD1", atom::AA_CD1, atom::C},
		{" CD2", atom::AA_CD2, atom::C},
		{" CE ", atom::AA_CE, atom::C},
		{" CE1", atom::AA_CE1, atom::C},
		{" CE2", atom::AA_CE2, atom::C},
		{" CE3", atom::AA_CE3, atom::C},
		{" CG ", atom::AA_CG, atom::C},
		{" CG1", atom::AA_CG1, atom::C},
		{" CG2", atom::AA_CG2, atom::C},
		{" CH2", atom::AA_CH2, atom::C},
		{" CZ ", atom::AA_CZ, atom::C},
		{" CZ2", atom::AA_CZ2, atom::C},
		{" CZ3", atom::AA_CZ3, atom::C},
		{" H  ", atom::AA_H, atom::H},
		{" H1 ", atom::AA_H1, atom::H},
		{" H2 ", atom::AA_H2, atom::H},
		{" H3 ", atom::AA_H3, atom::H},
		{" HA ", atom::AA_HA, atom::H},
		{" HA2", atom::AA_HA2, atom::H},
		{" HA3", atom::AA_HA3, atom::H},
		{" HB ", atom::AA_HB, atom::H},
		{" HB1", atom::AA_HB1, atom::H},
		{" HB2", atom::AA_HB2, atom::H},
		{" HB3", atom::AA_HB3, atom::H},
		{" HD1", atom::AA_HD1, atom::H},
		{"HD11", atom::AA_HD11, atom::H},
		{"HD12", atom::AA_HD12, atom::H},
		{"HD13", atom::AA_HD13, atom::H},
		{" HD2", atom::AA_HD2, atom::H},
		{"HD21", atom::AA_HD21, atom::H},
		{"HD22", atom::AA_HD22, atom::H},
		{"HD23", atom::AA_HD23, atom::H},
		{" HD3", atom::AA_HD3, atom::H},
		{" HE ", atom::AA_HE, atom::H},
		{" HE1", atom::AA_HE1, atom::H},
		{" HE2", atom::AA_HE2, atom::H},
		{"HE21", atom::AA_HE21, atom::H},
		{"HE22", atom::AA_HE22, atom::H},
		{" HE3", atom::AA_HE3, atom::H},
		{" HG ", atom::AA_HG, atom::H},
		{" HG1", atom::AA_HG1, atom::H},
		{"HG11", atom::AA_HG11, atom::H},
		{"HG12", atom::AA_HG12, atom::H},
		{"HG13", atom::AA_HG13, atom::H},
		{" HG2", atom::AA_HG2, atom::H},
		{"HG21", atom::AA_HG21, atom::H},
		{"HG22", atom::AA_HG22, atom::H},
		{"HG23", atom::AA_HG23, atom::H},
		{" HG3", atom::AA_HG3, atom::H},
		{" HH ", atom::AA_HH, atom::H},
		{"HH11", atom::AA_HH11, atom::H},
		{"HH12", atom::AA_HH12, atom::H},
		{" HH2", atom::AA_HH2, atom::H},
		{"HH21", atom::AA_HH21, atom::H},
		{"HH22", atom::AA_HH22, atom::H},
		{" HZ ", atom::AA_HZ, atom::H},
		{" HZ1", atom::AA_HZ1, atom::H},
		{" HZ2", atom::AA_HZ2, atom::H},
		{" HZ3", atom::AA_HZ3, atom::H},
		{" N  ", atom::AA_N, atom::N},
		{" ND1", atom::AA_ND1, atom::N},
		{" ND2", atom::AA_ND2, atom::N},
		{" NE ", atom::AA_NE, atom::N},
		{" NE1", atom::AA_NE1, atom::N},
		{" NE2", atom::AA_NE2, atom::N},
		{" NH1", atom::AA_NH1, atom::N},
		{" NH2", atom::AA_NH2, atom::N},
		{" NZ ", atom::AA_NZ, atom::N},
		{" O  ", atom::AA_O, atom::O},
		{" OD1", atom::AA_OD1, atom::O},
		{" OD2", atom::AA_OD2, atom::O},
		{" OE1", atom::AA_OE1, atom::O},
		{" OE2", atom::AA_OE2, atom::O},
		{" OG ", atom::AA_OG, atom::O},
		{" OG1", atom::AA_OG1, atom::O},
		{" OH ", atom::AA_OH, atom::O},
		{" OXT", atom::AA_OXT, atom::O},
		{" SD ", atom::AA_SD, atom::S},
		{" SG ", atom::AA_SG, atom::S},
		{" O  ", atom::AA_CEA, atom::O},
		{"UNKN", atom::N_INV, atom::T_INV}
	};

	static atom_data atom_name_nu_data[] = {
		// Nucleotides
		{" C1'", atom::NU_C1p, atom::C},{" C1*", atom::NU_C1p, atom::C},
		{" C2 ", atom::NU_C2, atom::C},
		{" C2'", atom::NU_C2p, atom::C},{" C2*", atom::NU_C2p, atom::C},
		{" C3'", atom::NU_C3p, atom::C},{" C3*", atom::NU_C3p, atom::C},
		{" C4 ", atom::NU_C4, atom::C},
		{" C4'", atom::NU_C4p, atom::C},{" C4*", atom::NU_C4p, atom::C},
		{" C5 ", atom::NU_C5, atom::C},
		{" C5'", atom::NU_C5p, atom::C},{" C5*", atom::NU_C5p, atom::C},
		{" C6 ", atom::NU_C6, atom::C},
		{" C7 ", atom::NU_C7, atom::C},
		{" C8 ", atom::NU_C8, atom::C},
		{" H1 ", atom::NU_H1, atom::H},
		{" H1'", atom::NU_H1p, atom::H},
		{" H2 ", atom::NU_H2, atom::H},
		{" H21", atom::NU_H21, atom::H},
		{" H22", atom::NU_H22, atom::H},
		{" H2'", atom::NU_H2p, atom::H},
		{"H2''", atom::NU_H2pp, atom::H},
		{" H3 ", atom::NU_H3, atom::H},
		{" H3'", atom::NU_H3p, atom::H},
		{" H41", atom::NU_H41, atom::H},
		{" H42", atom::NU_H42, atom::H},
		{" H4'", atom::NU_H4p, atom::H},
		{" H5 ", atom::NU_H5, atom::H},
		{" H5'", atom::NU_H5p, atom::H},
		{"H5''", atom::NU_H5pp, atom::H},
		{" H6 ", atom::NU_H6, atom::H},
		{" H61", atom::NU_H61, atom::H},
		{" H62", atom::NU_H62, atom::H},
		{" H71", atom::NU_H71, atom::H},
		{" H72", atom::NU_H72, atom::H},
		{" H73", atom::NU_H73, atom::H},
		{" H8 ", atom::NU_H8, atom::H},
		{"HO2'", atom::NU_HO2p, atom::H},
		{"HO3'", atom::NU_HO3p, atom::H},
		{"HO5'", atom::NU_HO5p, atom::H},
		{" N1 ", atom::NU_N1, atom::N},
		{" N2 ", atom::NU_N2, atom::N},
		{" N3 ", atom::NU_N3, atom::N},
		{" N4 ", atom::NU_N4, atom::N},
		{" N6 ", atom::NU_N6, atom::N},
		{" N7 ", atom::NU_N7, atom::N},
		{" N9 ", atom::NU_N9, atom::N},
		{" O2 ", atom::NU_O2, atom::O},
		{" O2'", atom::NU_O2p, atom::O},{" O2*", atom::NU_O2p, atom::O},
		{" O3'", atom::NU_O3p, atom::O},{" O3*", atom::NU_O3p, atom::O},
		{" O4 ", atom::NU_O4, atom::O},
		{" O4'", atom::NU_O4p, atom::O},{" O4*", atom::NU_O4p, atom::O},
		{" O5'", atom::NU_O5p, atom::O},{" O5*", atom::NU_O5p, atom::O},
		{" O6 ", atom::NU_O6, atom::O},
		{" OP1", atom::NU_OP1, atom::O},{" O1P", atom::NU_OP1, atom::O},
		{" OP2", atom::NU_OP2, atom::O},{" O2P", atom::NU_OP2, atom::O},
		{" OP3", atom::NU_OP3, atom::O},{" O3P", atom::NU_OP3, atom::O},
		{" O1A", atom::NU_O1A, atom::O},
		{" O2A", atom::NU_O2A, atom::O},
		{" O3A", atom::NU_O3A, atom::O},
		{" O1B", atom::NU_O1B, atom::O},
		{" O2B", atom::NU_O2B, atom::O},
		{" O3B", atom::NU_O3B, atom::O},
		{" O1G", atom::NU_O1G, atom::O},
		{" O2G", atom::NU_O2G, atom::O},
		{" O3G", atom::NU_O3G, atom::O},
		{" P  ", atom::NU_P, atom::P},
		{" PA ", atom::NU_PA, atom::P},
		{" PB ", atom::NU_PB, atom::P},
		{" PG ", atom::NU_PG, atom::P},
		{" O  ", atom::NU_CEN, atom::O},
		{"UNKN", atom::N_INV, atom::T_INV}
	};

	static atom_data atom_name_wa_data[] = {
		//Water
		{" O  ", atom::WA_O, atom::O},
		{"UNKN", atom::N_INV, atom::T_INV}
	};

	static atom_data atom_name_me_data[] = {
		{"MG  ", atom::ME_MG, atom::MG},
		{"ZN  ", atom::ME_ZN, atom::ZN},
		{"CA  ", atom::ME_CA, atom::CA},
		{"NA  ", atom::ME_NA, atom::NA},
		{"NI  ", atom::ME_NI, atom::NI},
		{"CL  ", atom::ME_CL, atom::CL},
		{"FE  ", atom::ME_FE, atom::FE},
		{"UNKN", atom::N_INV, atom::T_INV}
	};

};
