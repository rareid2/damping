IDIR =include
EIGEN_DIR=lib/eigen/
# Point to your Matlab install -- (specific to your machine)
MATLAB_DIR=/Applications/MATLAB_R2019a.app/extern/include
CC=c++
# CC=gcc

CFLAGS=-I$(IDIR) -I$(EIGEN_DIR) -I$(MATLAB_DIR) -L /usr/local/Cellar/gcc/8.3.0/lib/gcc/8 -L /usr/local/lib/ -lc++
# MATLAB_FLAGS = -L /usr/local/MATLAB/R2012a/bin/glnxa64 -leng -lm -lmx -lmex -lmat -lut -Wl,-rpath=/usr/local/MATLAB/R2012a/bin/glnxa64
MATLAB_FLAGS = -L /Applications/MATLAB_R2019a.app/bin/maci64 -leng -lm -lmx -lmex -lmat -lut -Wl -rpath /Applications/MATLAB_R2019a.app/bin/maci64

# g95
FORTRAN_FLAGS = -pg -Wall -fstatic -ffixed-line-length-132 -ffree-line-length-huge

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

# # The plasmapause location scripts, lifted from GCPM v2.4
# _FORTRAN_OBJ = pp_profile_d.mod types.mod util.mod constants.mod
# FORTRAN_OBJ = $(patsubst %,$(ODIR)/%,$(_FORTRAN_OBJ))

XFORM = lib/xform_double
# Rules for making individual objects
# $(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
# 	$(CC) -c -o $@ $< $(CFLAGS) -I$(EIGEN_DIR) -L$(LDIR) 
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I$(EIGEN_DIR)

# # Fortran modules
# $(ODIR)/%.mod: $(SRC_DIR)/%.f95 $(DEPS)
# 	g95 ${FORTRAN_FLAGS} -c -o $@ $<

# Rule to link everything together + generate executable
damping: $(OBJ) libxformd.a $(ODIR)/gauss_legendre.o
	# $(MAKE) -C $(XFORM)
	$(CC) $(CFLAGS) $(OBJ) $(ODIR)/gauss_legendre.o -L $(LDIR) -lxformd -lgfortran $(MATLAB_FLAGS) -o $(BDIR)/$@
	# Move the binaries up to the main directory
	cp bin/damping ../bin/damping
	cp data/crres_clean.mat ../bin/crres_clean.mat

dump: $(ODIR)/dump_psd_models.o $(ODIR)/kp_to_pp.o $(ODIR)/polyfit.o $(ODIR)/psd_model.o
	$(CC) $(CFLAGS) $(ODIR)/dump_psd_models.o $(ODIR)/kp_to_pp.o $(ODIR)/polyfit.o $(ODIR)/psd_model.o $(ODIR)/gauss_legendre.o -L $(LDIR) -lxformd -lgfortran $(MATLAB_FLAGS) -o $(BDIR)/$@

$(ODIR)/gauss_legendre.o: $(SRC_DIR)/gauss_legendre.c $(IDIR)/gauss_legendre.h

	gcc -c -o $@ $< $(CFLAGS)

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

