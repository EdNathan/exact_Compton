######################################################
#Makefile for 
# drive_SRF
#
#Variables:
# fc= Fortran compiler
# flags = general flags for compilation
# libs = libraries used during linking
# name = root name of the code
# code = Main fortran code
# exec = Binary executable (final output)
# objects = list of subroutines and functions used
#
# Use 'make' to compile everything
# Use 'make clean' to erase all the objects
# and 'make cleanout' to erase OUTPUT files.
#
#fc=ifort
#flags=-O3 -parallel -no-prec-div -shared-intel -xHost -fp-model source -mcmodel=large
#
# This worked on yorp, chandra, etc (UMD)
#flags=-O3 -parallel -no-prec-div -static-intel -xHost -fp-model source -mcmodel=large
#libs=$(HEADAS)/lib/libcfitsio_3.27.so
#
#flags= -O3 -fp-model source -fopenmp
#flags= -O3 -check bounds -fp-model source
#


flags=-O3
#flags=-O3 -fallow-argument-mismatch -fopenmp
#flags=-O3 -fallow-argument-mismatch
#flags=-Og -g -fcheck=all -Wextra -fimplicit-none -fbacktrace -fallow-argument-mismatch 

ifeq ($(c), intel)
 fc=ifx
 # For debugging:
 # flags=-O0 -g -traceback # -check all -debug all
 flags += -mcmodel=large -shared-intel
else
 fc=gfortran
endif



p ?= f

# If Make is called with p=t, add the -fopenmp
# This allows the compiler to understand the OMP library for parallisation
# Note (for future):  -fopenmp is applicable for gfortran.  Intel Compilers use -qopenmp
ifeq ($(p), t)
 ifeq ($(c), intel)
  flags += -qopenmp
 else
  flags += -fopenmp
 endif
endif

libs = -lcfitsio

#flags=-march=native -ffast-math -funroll-loops -O3 -finline-limit=600
#flags= -O3 -fbounds-check
#
#libs= -lcfitsio
name=drive_SRF
code=$(name).f
exec=$(name).x
mysrc=my-routines
mymodsrc=my-modules

myobjts= $(mymodsrc)/constants.o                 \
         $(mymodsrc)/fits_writing.o              \
         $(mysrc)/srf_nonlimit.o                 \
         $(mysrc)/bk2.o                          \
         $(mysrc)/crsexact.o                     \
         $(mysrc)/enegrd.o                       \
         $(mysrc)/gaulegf.o           	         \
         $(mysrc)/probab.o           	         \
         $(mysrc)/scattxs.o                      \
         $(mysrc)/super_Compton_RF_fits.o        \
         $(mysrc)/super_Compton_RF_fits_angle.o  \


# Compile xstar with all subroutines
$(exec):  $(myobjts)
	$(fc) $(flags) $(code) $(myobjts) $(libs) -o $(exec)

# Compile and create objects (XILLVER)
$(myobjts): %.o: %.f
	$(fc) $(flags) -c $< -o $@

# Clean all objects
clean:
	rm $(mysrc)/*.o
	rm $(mymodsrc)/*.o
	rm *.mod
	rm $(exec)

all:
	rm $(exec)
	make $(exec)

# Tue Dec 18 14:13:06 EST 2007
# Javier Garcia
# Modified Dec 28 T. Kallman
# Modified Dec 01 (2022) E. Nathan
