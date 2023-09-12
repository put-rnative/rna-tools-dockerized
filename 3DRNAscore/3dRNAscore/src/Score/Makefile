CXX      := g++
CXXFLAGS  := -std=c++11
LD      := ld
LDFLAGS := 
SOURCE  := $(wildcard *.cpp)
OBJS    := $(patsubst %.cpp, %.o, $(SOURCE))

all:$(OBJS)

%.o:%.cpp
	@echo Compiling $< ...
	$(CXX) -c $(CXXFLAGS) $< -o $*.o
 
clean:
	@rm -f *.o 

