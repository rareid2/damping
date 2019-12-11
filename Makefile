IDIR =include
EIGEN_DIR=lib/eigen/
# Point to your Matlab install -- (specific to your machine)
MATLAB_DIR=/Applications/MATLAB_R2019a.app/extern/include
CC=c++


CFLAGS=-I$(IDIR) -I$(EIGEN_DIR) -I$(MATLAB_DIR) -L /usr/local/Cellar/gcc/8.3.0/lib/gcc/8 -L /usr/local/lib/ -lc++

# Nansen flags
# MATLAB_FLAGS = -L /usr/local/MATLAB/R2012a/bin/glnxa64 -leng -lm -lmx -lmex -lmat -lut -Wl,-rpath=/usr/local/MATLAB/R2012a/bin/glnxa64
# OSX flags
MATLAB_FLAGS = -L /Applications/MATLAB_R2019a.app/bin/maci64 -leng -lm -lmx -lmex -lmat -lut -Wl -rpath /Applications/MATLAB_R2019a.app/bin/maci64

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

# Objects to build for dump_psd_models: (both damping_main and dump_psd_models have 'main' functions)
_DUMP_OBJ = dump_psd_models.o kp_to_pp.o polyfit.o psd_model.o
DUMP_OBJ = $(patsubst %,$(ODIR)/%,$(_DUMP_OBJ))

XFORM = lib/xform_double
# Rules for making individual objects
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I$(EIGEN_DIR)

# Rule to link everything together + generate executable
# (This is the main executable, "damping")
damping: $(OBJ) libxformd.a $(ODIR)/gauss_legendre.o
	$(CC) $(CFLAGS) $(OBJ) $(ODIR)/gauss_legendre.o -L $(LDIR) -lxformd -lgfortran $(MATLAB_FLAGS) -o $(BDIR)/$@
	# Move the binaries up to the main directory
	cp bin/damping ../bin/damping
	cp data/crres_clean.mat ../bin/crres_clean.mat

# Link and build "dump_psd_models" executable 
dump_psd_models: $(DUMP_OBJ) libxformd.a
	$(CC) $(CFLAGS) $(DUMP_OBJ) -L $(LDIR) -lxformd -lgfortran $(MATLAB_FLAGS) -o $(BDIR)/$@ 

# Need a separate rule for gauss_legendre (I forget why)
$(ODIR)/gauss_legendre.o: $(SRC_DIR)/gauss_legendre.c $(IDIR)/gauss_legendre.h
	gcc -c -o $@ $< $(CFLAGS)

# Make the coordinate transformation library
libxformd.a:
	$(MAKE) -C $(XFORM)

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

