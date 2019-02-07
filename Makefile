##########################
# Macros: name = data
##########################

CXX=g++
CXXFLAGS= -Wall -Wno-deprecated  -g -I$(HOME) -O0 -Wextra -pedantic -g3 -L -lboost_thread -fopenmp

CXXFLAGS_PROFILE=-Wall -Wno-deprecated  -pg -I$(HOME) -O4 -Wextra -pedantic -L -lboost_thread -fopenmp

CXXFLAGS = $(CXXFLAGS_PROFILE)
CCFLAGS  = $(CXXFLAGS_PROFILE)

LIBS= -lm
LDFLAGS = -pg

#########################
#Rules: target: source
#		command
#########################

#.SUFFIXES: .cc .c .C

#.C.o:
#	$(CC) 	-c 	$(CFLAGS) 	-o 	$@ 	$<

all: 	hYrnpt hYrnpt_par hYrnpt_msd hYrnpt_parmsd

hYrnpt:	ran2.o tools.o inout.o setup.o hYrnpt.o lineDist.o 
	$(CXX) $(LDFLAGS) $(LIBS) ran2.o lineDist.o tools.o inout.o setup.o hYrnpt.o -o hYr.exe -g -lstdc++ -L -lboost_thread -fopenmp

hYrnpt_sh:	lineDist.o ran2.o tools.o inout.o setup.o hYrnpt_shape.o 
	$(CXX) $(LDFLAGS) $(LIBS) ran2.o lineDist.o tools.o inout.o setup.o hYrnpt_shape.o -o hYr_sh.exe -g -lstdc++ -L -lboost_thread -fopenmp

hYrnpt_par: lineDist.o ran2.o tools.o inout.o setup.o hYrnpt_par.o
	$(CXX) $(LDFLAGS) $(LIBS) ran2.o lineDist.o tools.o inout.o setup.o hYrnpt_par.o -o hYr_par.exe -g -lstdc++ -L -lboost_thread -fopenmp

hYrnpt_msd:	lineDist.o ran2.o tools.o inout.o setup.o hYrnpt_msd.o 
	$(CXX) $(LDFLAGS) $(LIBS) ran2.o lineDist.o tools.o inout.o setup.o hYrnpt_msd.o -o hYr_msd.exe -g -lstdc++ -L -lboost_thread -fopenmp

hYrnpt_parmsd: lineDist.o ran2.o tools.o inout.o setup.o hYrnpt_parmsd.o
	$(CXX) $(LDFLAGS) $(LIBS) ran2.o lineDist.o tools.o inout.o setup.o hYrnpt_parmsd.o -o hYr_parmsd.exe -g -lstdc++ -L -lboost_thread -fopenmp

.PHONY: clean
clean: 
	rm -f *~ \#*
	rm -f hYrnpt.o lineDist.o ran2.o tools.o inout.o setup.o hYrnpt_shape.o hYrnpt_par.o hYrnpt_msd.o hYrnpt_parmsd.o
