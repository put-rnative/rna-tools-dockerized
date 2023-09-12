CXX        := g++
CXXFLAGS   := -std=c++11

TOOLS_H  := $(wildcard Tools/*.h)
TOOLS_CPP  := $(wildcard Tools/*.cpp)
TOOLS_OBJ  := $(patsubst %.cpp,%.o,$(TOOLS_CPP))
PDB_H    := $(wildcard Pdb/*.h)
PDB_CPP    := $(wildcard Pdb/*.cpp)
PDB_OBJ    := $(patsubst %.cpp,%.o,$(PDB_CPP))
SCORE_H  := $(wildcard Score/*.h)
SCORE_CPP  := $(wildcard Score/*.cpp)
SCORE_OBJ  := $(patsubst %.cpp,%.o,$(SCORE_CPP))

all:3dRNAscore/3dRNAscore-1.0

3dRNAscore/3dRNAscore-1.0:lib3dRNA.a 3dRNAscore/3dRNAscore.o
	$(CXX) $(CXXFLAGS) -o 3dRNAscore/3dRNAscore-1.0 3dRNAscore/3dRNAscore.o lib3dRNA.a

3dRNAscore/3dRNAscore.o:3dRNAscore/3dRNAscore.cpp $(TOOLS_OBJ) $(PDB_OBJ) $(SCORE_OBJ)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

lib3dRNA.a:$(TOOLS_OBJ) $(PDB_OBJ) $(SCORE_OBJ)
	ar -rc lib3dRNA.a $(TOOLS_OBJ) $(PDB_OBJ) $(SCORE_OBJ)

Tools/%.o:Tools/%.cpp Tools/%.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

Pdb/%.o:Pdb/%.cpp Pdb/%.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

Score/%.o:Score/%.cpp Score/%.h
	$(CXX) -c $(CXXFLAGS) -o $@ $<

clean:
	@rm -f $(TOOLS_OBJ) $(PDB_OBJ) $(SCORE_OBJ) lib3dRNA.a 3dRNAscore/3dRNAscore.o 3dRNAscore/3dRNAscore-1.0

