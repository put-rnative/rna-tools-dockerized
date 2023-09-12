CXX      := g++
CXXFLAGS  :=  -std=c++11
LD      := ld
LDFLAGS := 
SOURCE  := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp, %.o, $(SOURCE))

all:3dRNAscore-1.0

3dRNAscore-1.0:$(OBJS) ../lib3dRNA.a
	$(CXX) $(CXXFLAGS) $(OBJS) ../lib3dRNA.a -o 3dRNAscore-1.0

%.o:%.cpp
	@echo Compiling $< ...
	$(CXX) -c $(CXXFLAGS) $< -o $*.o
 
clean:
	@rm -f *.o 3dRNAscore-1.0

