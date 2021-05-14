# Author SHI Xin <shixin@ihep.ac.cn>  
# Created [2021-03-31 Sun 14:38] 

MYPWD=${PWD}

BIN=./bin
SRC=./src

CC = g++
GCCFLAGS  = -Wall -g 

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS))

#FLAGS=$(GCCFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -I $(MYPWD)/include/  -I /usr/local/root/include/ -I /usr/include/eigen3/ -L/usr/local/root/lib -L $(MYPWD)/lib -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib -ldolfin -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -lTreePlayer -lTreeViewer -lHistPainter  -lFFTW -lFITSIO -lGX11TTF  -lMinuit2 -lMathMore -lCling -lRooFit -lRooFitCore -lMatrix -lboost_system
#yuhang 2021.5.14 yuhang
FLAGS=$(GCCFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -lHistPainter 

# 2021.5.4 jia ning 
#FLAGS=$(GCCFLAGS) $(ROOTCFLAGS) $(ROOTLIBS) -I $(MYPWD)/include/  -I /usr/local/root/include/ -I /usr/include/eigen3/ -L/usr/local/root/lib -L $(MYPWD)/lib -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib -L/usr/local/lib -ldolfin -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -lTreePlayer -lTreeViewer -lHistPainter  -lFFTW -lFITSIO -lGX11TTF  -lMinuit2 -lMathMore -lCling -lRooFit -lRooFitCore -lMatrix -lboost_system


PROG=raser
LIST=$(addprefix $(BIN)/, $(PROG))


all: $(LIST)
	@echo "Build successful."

$(LIST): | $(BIN)

$(BIN): 
	mkdir -p $(BIN)

# $(BIN)/raser: $(SRC)/raser.cc
# 	$(CC) $< $(FLAGS) -o $@
$(BIN)/raser: $(SRC)/raser.cc
	$(CC) $< $(FLAGS) -o $@
merge: 
	git remote update && git merge dt-np/main 

clean:
	rm -f $(BIN)/raser 
	