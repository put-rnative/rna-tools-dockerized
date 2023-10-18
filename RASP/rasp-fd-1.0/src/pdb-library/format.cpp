#include <string>
#include <iostream>
#include <set>
#include "format.h"

namespace tnpdb {

	//////////////////////////////////////////////////////////////////////////////////////////////
	record record_line(const char *line) {
		if (line[0] == '\0') return R_OTHER;
		
		std::string rec(line,0,6);
		
		if (rec == "DBREF ") return R_DBREF;
		else if (rec == "SEQRES") return R_SEQRES;
		else if (rec == "ATOM  ") return R_ATOM;
		else if (rec == "HETATM") return R_HETATM;
		else if (rec == "MASTER") return R_MASTER;
		else if (rec == "ENDMDL") return R_ENDMDL;
		else if (rec == "END   " || rec == "END  " || rec == "END " || rec == "END") return R_END;
		else if (rec == "HEADER" || rec == "TITLE " || rec == "COMPND" || rec == "SOURCE" ||
			rec == "KEYWDS" || rec == "EXPDTA" || rec == "AUTHOR" || rec == "REVDAT" ||
			rec == "JRNL  " || rec == "REMARK" || rec == "HELIX " || rec == "SHEET " ||
			rec == "SITE  " || rec == "CRYST1" || rec == "ORIGX1" || rec == "ORIGX2" ||
			rec == "ORIGX3" || rec == "SCALE1" || rec == "SCALE2" || rec == "SCALE3" ||
			rec == "SEQADV" || rec == "TURN  " || rec == "FORMUL" || rec == "HETNAM" ||
			rec == "SSBOND" || rec == "MODRES" || rec == "HET   " || rec == "CISPEP") return R_HEADER;
		else if (rec == "CONECT") return R_CONECT;
		else if (rec == "TER   " || rec == "TER") return R_TER;
		else if (rec == "MODEL ") return R_MODEL;
		else {
			std::cerr << "\"" + rec + "\" is not a known line type." << std::endl;
			return R_OTHER;
		}
	}

};
