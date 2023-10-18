namespace tnpdb {
	
	enum _record {R_HEADER, R_DBREF, R_SEQRES, R_ATOM, R_HETATM, R_MASTER, R_ENDMDL, R_OTHER, R_TER, R_MODEL, R_CONECT, R_END };
	typedef _record record;

	const char
		atom_line_iformat[] = "ATOM  %5d%*1c%4c%1c%3c%*c%1c%4d%1c%*1c%*1c%*1c%8f%8f%8f%6f%6f%*1c%*1c%*1c%*1c%*1c%*1c%4c%2c%2c",
		atom_line_oformat[] = "ATOM  %5d %4s%1c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s",
	    ter_line_oformat[] = "TER   %5d      %3s %c%4d%c",
		hetatm_line_iformat[] = "HETATM%5d%*1c%4c%1c%3c%*c%1c%4d%1c%*1c%*1c%*1c%8f%8f%8f%6f%6f%*1c%*1c%*1c%*1c%*1c%*1c%4c%2c%2c",
		hetatm_line_oformat[] = "HETATM%5d %4s%1c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s";
	
	record record_line(const char *);

};
