include make.inc

LIB        = libprolate
FC         = gfortran
ifeq ($(FC),gfortran)
    FFLAGS = -O3 -Wall
else ifeq ($(FC),ifort)
    FFLAGS = -O3 -warn all
else
    FFLAGS = -O3
endif

ARCH       = ar
ARCHFLAGS  = cr
RANLIB     = ranlib

SRC        = src
EXAMPLE    = example

vpath %.f90 $(SRC)

OBJ        = prolate.o
EXEC       = int2

LDLIBS     = -lblas -llapack -lprint

.PHONY: all driver

all: $(LIB).a

$(LIB).a: $(OBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(OBJ)
	$(RANLIB) $@

driver: $(LIB).a
	$(FC) $(FFLAGS) -I$(INCDIR) -L$(LIBDIR) -L. -o $(EXEC) $(EXAMPLE)/prolate_driver.f90 $(LDLIBS) -lprolate

clean:
	rm -f *.o

cleanlib:
	rm -f *.mod
	rm -f $(LIB).a

cleandriver:
	rm -f $(EXEC)

cleanmat:
	rm -f matlab/*.mod
	rm -f matlab/*.mex*

cleanall: clean cleanlib cleandriver cleanmat

%.o: %.f90
	$(FC) $(FFLAGS) -I$(INCDIR) -c $< -o $@

help:
	@echo "Please use \"make <target>\" where <target> is one of"
	@echo "  all         to make $(LIB).a library"
	@echo "  driver      to make the example driver program"
	@echo "  clean       to remove all compiled objects"
	@echo "  cleanlib    to remove all libraries"
	@echo "  cleandriver to remove example driver program"
	@echo "  cleanmat    to remove the matlab mex files"
	@echo "  cleanall    to clean everything"
	@echo "  help        to display this help message"
