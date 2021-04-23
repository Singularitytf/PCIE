# --------------------------------------------------------------
# GNUmakefile for ROOT.                      Liuhm, June 1, 2003
# Update for multi-file code system	 C.-h. Fang Apr 18, 2021
# --------------------------------------------------------------

#CPPFLAGS += -I/usr/local/include/root        #Change!
CPPFLAGS += $(shell root-config --cflags)    #Change!
LDLIBS += -O3
LDLIBS += $(shell root-config --glibs)

#-lRooFit
#LDFLAGS:= -L/usr/local/include/root         #Change!
LDFLAGS:= -L$(shell root-config --libdir)         #Change!
SRCDIR += ./src

SRCs		= ./src/debug.cc ./src/Belt.cc ./src/Poisson.cc  #Change!
OBJs		= debug.o Belt.o Poisson.o
# Short for Poisson Confident Interval Estimator.
PROGRAM		= PCIE

cpl:	$(OBJs)
link:	$(PROGRAM)
all:	$(OBJs) $(PROGRAM)

$(PROGRAM):	$(OBJs)
	@echo "Linking $(PROGRAM) ..."
	@echo "g++ $(OBJs) $(CPPFLAGS) $(LDLIBS) -o $(PROGRAM)"
	@g++ $(OBJs) $(CPPFLAGS) $(LDLIBS) -o $(PROGRAM)
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


