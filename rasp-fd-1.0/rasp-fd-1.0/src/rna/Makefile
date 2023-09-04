CC = g++
PDBLIB = ../pdb-library
TOOLIB = ../tool-library

all: rasp_fd rasp_profile_fd

rasp_fd: rasp_fd.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o
	${CC} ${CFLAG} rasp_fd.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o -o rasp_fd -lm

rasp_profile_fd: rasp_profile_fd.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o
	${CC} ${CFLAG} rasp_profile_fd.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o -o rasp_profile_fd -lm

rna_frequency_matrix: rna_frequency_matrix.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o
	${CC} ${CFLAG} rna_frequency_matrix.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o -o rna_frequency_matrix -lm

rna_potential_matrix: rna_potential_matrix.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o
	${CC} ${CFLAG} rna_potential_matrix.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o -o rna_potential_matrix -lm

rna_potential_eval: rna_potential_eval.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o
	${CC} ${CFLAG} rna_potential_eval.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o -o rna_potential_eval -lm

rna_potential_profile: rna_potential_profile.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o
	${CC} ${CFLAG} rna_potential_profile.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o -o rna_potential_profile -lm

rna_potential_plots: rna_potential_plots.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o
	${CC} ${CFLAG} rna_potential_plots.o ${PDBLIB}/pdb.o ${PDBLIB}/model.o ${PDBLIB}/chain.o ${PDBLIB}/residue.o ${PDBLIB}/atom.o ${PDBLIB}/format.o ${TOOLIB}/string.o ${TOOLIB}/matrix.o rna_potential.o -o rna_potential_plots -lm

rna_frequency_matrix.o: rna_frequency_matrix.cpp ${PDBLIB}/atom.h
	${CC} ${CFLAG} -c rna_frequency_matrix.cpp -o rna_frequency_matrix.o
	
rna_potential_matrix.o: rna_potential_matrix.cpp ${PDBLIB}/atom.h
	${CC} ${CFLAG} -c rna_potential_matrix.cpp -o rna_potential_matrix.o

rna_potential_eval.o: rna_potential_eval.cpp ${PDBLIB}/atom.h
	${CC} ${CFLAG} -c rna_potential_eval.cpp -o rna_potential_eval.o

rna_potential_profile.o: rna_potential_profile.cpp ${PDBLIB}/atom.h ${TOOLIB}/gradient.h
	${CC} ${CFLAG} -c rna_potential_profile.cpp -o rna_potential_profile.o
		
rna_potential_plots.o: rna_potential_plots.cpp ${PDBLIB}/atom.h
	${CC} ${CFLAG} -c rna_potential_plots.cpp -o rna_potential_plots.o

rna_potential.o: rna_potential.cpp rna_potential.h ${TOOLIB}/string.h
	${CC} ${CFLAG} -c rna_potential.cpp -o rna_potential.o
	
rasp_fd.o: rasp_fd.cpp ${PDBLIB}/atom.h
	${CC} ${CFLAG} -c rasp_fd.cpp -o rasp_fd.o
	
rasp_profile_fd.o: rasp_profile_fd.cpp ${PDBLIB}/atom.h ${TOOLIB}/gradient.h
	${CC} ${CFLAG} -c rasp_profile_fd.cpp -o rasp_profile_fd.o

clear:
	\rm *.o

