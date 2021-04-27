# --------------------------------------------------------------
# GNUmakefile for ROOT.                      Liuhm, June 1, 2003
# Update for multi-file code system	 C.-h. Fang Apr 18, 2021
# --------------------------------------------------------------

#CPPFLAGS += -I/usr/local/include/root        #Change!
CPPFLAGS += $(shell root-config --cflags)    #Change!
CPPFLAGS += $(shell gsl-config --cflags)
LDLIBS += -O3
LDLIBS += $(shell root-config --glibs)
LDLIBS += $(shell gsl-config --libs)

#-lRooFit
#LDFLAGS:= -L/usr/local/include/root         #Change!
LDFLAGS:= -L$(shell root-config --libdir)         #Change!
SRCDIR += ./src

SRCs		= ./src/debug.cc ./src/Interval.cc ./src/Belt.cc ./src/Poisson.cc  ./src/cvlPoisson.cc #Change!
OBJs		= debug.o Interval.o Belt.o Poisson.o cvlPoisson.o
# Short for Poisson Confident Interval Estimator.
PROGRAM		= PCIE

cpl:	$(OBJs)
link:	$(PROGRAM)
all:	$(OBJs) $(PROGRAM)

$(PROGRAM):	$(OBJs)
	@echo "Linking $(PROGRAM) ..."
	@echo "g++ $(OBJs) $(CPPFLAGS) $(LDLIBS) -o $(PROGRAM)"
	@g++ $(OBJs) $(CPPFLAGS) $(LDLIBS) -lgsl -o $(PROGRAM)
	@echo "done"
$(OBJs):        $(SRCs)
	@echo "Compiling $(OBJs)...."
	@echo "g++ -c $(CPPFLAGS) $(SRCs)"
	@g++ -c $(CPPFLAGS) $(SRCs)


#clean: 	$(PROGRAM)
#	@echo "rm -f $(OBJ1) $(OBJ2) core"
#	@rm -f $(OBJ1) $(OBJ2) outputfilename.h.gch $(PROGRAM)


#edit : writecode.o 
#	g++ -o edit writecode.o
 
#writecode.o : writecode.cc 
#	g++ writecode.cc

clean :
	@echo "rm -f $(OBJs) *.h.gch $(PROGRAM)"
	@rm -f $(OBJ1) *.h.gch $(PROGRAM) *.o


