IDIR =include
EIGEN_DIR=lib/eigen-3.3.7/

# Point to your Matlab install -- (specific to your machine)
MATLAB_DIR=/Applications/MATLAB_R2019a.app/extern/include
MATLAB_BIN=/Applications/MATLAB_R2019a.app/bin/maci64

# Path to your libgfortran library, wherever it may reside. /usr/lib? /usr/local/lib? 
# (Mine is in Cellar since I installed it using Homebrew)
LGFORTRAN_PATH=/usr/local/Cellar/gcc/8.3.0/lib/gcc/8

# Which compiler to use -- c++, gcc, etc
CC=c++

CFLAGS=-I$(IDIR) -I$(EIGEN_DIR) -I$(MATLAB_DIR)

# Nansen flags
# MATLAB_FLAGS = -L /usr/local/MATLAB/R2012a/bin/glnxa64 -leng -lm -lmx -lmex -lmat -lut -Wl,-rpath=/usr/local/MATLAB/R2012a/bin/glnxa64
# OSX flags
MATLAB_FLAGS = -L $(MATLAB_BIN) -leng -lm -lmx -lmex -lmat -lut -Wl -rpath $(MATLAB_BIN)

# compiled module directory
ODIR =build

# Libraries
LDIR =lib
	
# output binary directory
BDIR =bin

# source files here
SRC_DIR=src

# Dependencies (header files)
_DEPS = damping.h consts.h coord_transforms.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# Objects to build
_OBJ = damping_main.o damping_ngo.o damping_foust.o math_utils.o \
	   kp_to_pp.o polyfit.o psd_model.o integrand.o wipp_fileutils.o \
	   bulge.o coord_transforms.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Objects to build for dump_psd_models: (we need separate lists, since
# both damping_main and dump_psd_models have 'main' functions)
_DUMP_OBJ = dump_psd_models.o kp_to_pp.o polyfit.o psd_model.o
DUMP_OBJ = $(patsubst %,$(ODIR)/%,$(_DUMP_OBJ))

# The Fortran coordinate transformation library
XFORM = lib/xform_double

# test -d $(ODIR) || mkdir $(ODIR)
# test -d $(BDIR) || mkdir $(BDIR)
# Rules for making individual objects
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I$(EIGEN_DIR)

# Rule to make all
all: $(ODIR) $(BDIR) damping dump_psd_models

# Rule to link everything together + generate executable
# (This is the main executable, "damping")
damping: $(OBJ) libxformd.a $(ODIR)/gauss_legendre.o
	test -d $(ODIR) || mkdir $(ODIR)
	test -d $(BDIR) || mkdir $(BDIR)
	$(CC) $(CFLAGS) $(OBJ) $(ODIR)/gauss_legendre.o -L $(LDIR) -L $(LGFORTRAN_PATH)\
								 -lxformd -lgfortran $(MATLAB_FLAGS) -o $(BDIR)/$@
	# Move the binaries up to the main directory
	# -- this assumes 'damping' is a subdirectory of the Stanford Raytracer
	cp bin/damping ../bin/damping
	cp data/crres_clean.mat ../bin/crres_clean.mat

# Link and build "dump_psd_models" executable 
dump_psd_models: $(DUMP_OBJ) libxformd.a
	$(CC) $(CFLAGS) $(DUMP_OBJ) -L $(LDIR) -L $(LGFORTRAN_PATH) -lxformd -lgfortran $(MATLAB_FLAGS) -o $(BDIR)/$@ 

# Need a separate rule for gauss_legendre (I forget why)
$(ODIR)/gauss_legendre.o: $(SRC_DIR)/gauss_legendre.c $(IDIR)/gauss_legendre.h
	gcc -c -o $@ $< $(CFLAGS)

# Make the coordinate transformation library
libxformd.a:
	$(MAKE) -C $(XFORM)
# 	cp $(XFORM)/libxformd.a $(LDIR)/libxformd.a

# Safety rule for any file named "clean"
.PHONY: clean

# Purge the build and bin directories
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	$(MAKE) -C $(XFORM) clean

# Everything except the transform library, since it hasn't changed
tidy:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	# $(MAKE) -C $(XFORM) clean

# Make the output directories, if they don't already exist
$(ODIR):
	test -d $(ODIR) || mkdir $(ODIR)
$(BDIR):
	test -d $(BDIR) || mkdir $(BDIR)

