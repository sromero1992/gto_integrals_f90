CC     = gcc
F90    = gfortran

LAPACK =Y
DEBUG  =Y

LIBS = 
ifeq ($(LAPACK),Y)
  LIBS = -llapack -lblas -lm
endif

ifeq ($(DEBUG),Y)
  # Debug flags
  #FFLAGS = -O3 -fbacktrace -fbounds-check -Wall -std=f2008 -fall-intrinsics -g -ffpe-trap=zero,overflow,underflow  #-funroll-loops 
  FFLAGS = -O0 -fbacktrace -fbounds-check -Wall -std=f2008 -fall-intrinsics -g -ffpe-trap=zero,overflow,underflow  #-funroll-loops 
else
  # Production flags
  FFLAGS = -O3 -fbacktrace 
endif

.SUFFIXES: .c .cc .ftn .f .o .f90 

# Compiling rule
%.o : %.f90
	$(F90) $(FFLAGS)  -c  $< 

# Object modules
modules = modules_diag.o modules_com.o basis_functions.o modules_g.o
#mainp   = diagon.o diagon_fortran.o  dens.o sort_eivec_val.o quick_sort.o read_info.o  okven_matbuild.o coul_build.o s_transform.o get_energy.o nuc_rep.o hf_scf.o basis_moments.o  getRmoment.o degsym.o  main_hf.o
mainp   = diagon.o dens.o sort_eivec_val.o quick_sort.o read_info.o  okven_matbuild.o coul_build.o s_transform.o get_energy.o nuc_rep.o hf_scf.o basis_moments.o  getRmoment.o degsym.o  main_hf.o
obj = $(modules) $(mainp)

bin = hf_scf.exe

# Default target
all : hf_scf

hf_scf : $(obj) 
	$(F90) $(FFLAGS) $(obj) -o  $(bin) $(LIBS) 

clean: 
	rm -f *.o *.mod *.txt
