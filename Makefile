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
fc=gfortran

#flags=-O3
#flags=-O3 -fallow-argument-mismatch -fopenmp
flags=-O3 -fallow-argument-mismatch
#flags=-Og -fcheck=all -Wextra -fimplicit-none -fbacktrace -fallow-argument-mismatch 

# If Make is called with the -k flag, add the -fopenmp
# This allows the compiler to understand the OMP library for parallisation
# Note (for future):  -fopenmp is applicable for gfortran.  Intel Compilers use -qopenmp
ifneq (,$(findstring k,$(MAKEFLAGS)))
flags += -fopenmp
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

myobjts= $(mysrc)/constants.o                    \
         $(mysrc)/bk2.o                          \
         $(mysrc)/crsexact.o                     \
         $(mysrc)/enegrd.o                       \
         $(mysrc)/gaulegf.o           	         \
         $(mysrc)/probab.o           	         \
         $(mysrc)/scattxs.o                      \
         $(mysrc)/write_fits.o                   \
         $(mysrc)/write_fits_angle.o             \
         $(mysrc)/super_Compton_RF.o             \
         $(mysrc)/super_Compton_RF_fits.o        \
         $(mysrc)/super_Compton_RF_fits_angle.o  \
         


# Compile xstar with all subroutines
$(exec): $(myobjts)
	$(fc) $(flags) $(code) $(myobjts) $(libs) -o $(exec)

# Compile and create objects (XILLVER)
$(myobjts): %.o: %.f
	$(fc) $(flags) -c $< -o $@

# Clean all objects
clean:
	rm $(mysrc)/*.o
	rm *.mod
	rm $(exec)

all:
	rm $(exec)
	make $(exec)

# Tue Dec 18 14:13:06 EST 2007
# Javier Garcia
# Modified Dec 28 T. Kallman
# Modified Dec 01 (2022) E. Nathan