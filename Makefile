# This is the makefile to compile the file main.f90 
SOURCE0 = main.f90 
SOURCE +=        			

MACRO+=-D__solve_topological_charge
#MACRO+=-D__solve_eigenstates

exec=main.out

OBJ0 = $(patsubst %.f90,%.o,$(SOURCE0))
OBJ = $(patsubst %.f90,%.o,$(SOURCE))

F90 = ifort 



FFLAGS_FPP =  
FFLAGS_FPP += -fpp 

LDFLAGS += -L/usr/local/intel/mkl/10.0.5.025/lib/em64t/ -static
FLIBS += -lmkl_lapack -lmkl_em64t -lmkl_solver -lguide -lpthread



#.SURFFIXES:

%.o : %.F90
	$(F90) $(FFLAGScompile) $(FFLAGS)  $(MACRO) $(FFLAGS_FPP) -c -o $@ $<
%.o : %.f90
	$(F90) $(FFLAGScompile) $(FFLAGS)  $(MACRO) $(FFLAGS_FPP) -c -o $@ $<

all:	$(exec)

$(OBJ) : $(OBJ0)

$(exec):	$(OBJ) $(OBJ0)
	$(F90) $(FFLAGSlink) $(LDFLAGS) -o $@ $^ $(FLIBS)


etags:TAGS
TAGS:	$(SOURCE)  $(SOURCE0)
	etags $^

clean:	
	rm -f *.od *.o  *_mod.f90 *.dat
	rm -f $(exec) nohup.out
  

distclean:
	rm -f *.od *.o *.mod *~ 
	rm -f $(exec)
	rm -f *.dat nohup.out PI*
