ROOTLIBS:= $(shell root-config --libs)

name=readWD

all: $(name)

clean:  
	rm *.o $(name)

$(name).o: $(name).cpp
	g++ -g3 -Wall -I $(ROOTSYS)/include/root -I $(ROOTSYS)/include -std=c++17 -o $(name).o -c $(name).cpp

$(name): $(name).o
	g++ -g3 $(name).o $(ROOTLIBS) -lncurses -o $(name) -Wall
