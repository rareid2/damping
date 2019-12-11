# damping
A tool to compute wave attenuation due to Landau damping and geometric spreading in the magnetosphere. 
A companion code to the Stanford 3D Ray Tracing package.

# Building
This software has been compiled on OSX Mojave (10.14.5), and Linux (Fedora 6), but should build on anything.

This software has the following requirements:
-- Eigen (http://eigen.tuxfamily.org). A zipped version is provided in ```lib/eigen-3.3.7.zip```.
Unzip it into the ```lib/``` directory.

-- xform_double: A set of Fortran coordinate transforms, originating with the Tsyganenko Geopack libraries. 
It's included in lib/xform_double.

-- Matlab: The code needs Matlab, and the associated C headers, to read a .mat data file. Edit the root Makefile to
point to your Matlab install, e.g.:
```
  MATLAB_DIR=/Applications/MATLAB_R2019a.app/extern/include
  MATLAB_BIN=/Applications/MATLAB_R2019a.app/bin/maci64
```

Depending on your machine, you may need to edit the Makefile to point the linker to your installed version of ```libgfortran```.
Mine shows up under ```LGFORTRAN_PATH=/usr/local/Cellar/gcc/8.3.0/lib/gcc/8```, since I installed gfortran as part of
gcc-8, using Homebrew.

With all the requirements satisfied, build the sofware by typing ```make``` from the root directory.

# Using
There are two executables generated: ```damping``` and ```dump_psd_models```.

## Damping
This is the main Landau damping code. It is run from the command line (```./damping```) with the following parameters 
passed as command-line arguments:
```
-------- Stanford Ray Tracer Landau Damping Code ---------- 
--inp_file, -i:     Input file to use  (default: input.ray)
--out_file, -o:     Output file to use (default: output.ray)
--mode, -m:         Mode -- 1 for modern implementation, 0 for the legacy
                    damping code from the 2d raytracer.
--AE, -a:           AE index - integer valued, 1, 2, 3.
--Kp, -k:           Kp: Real valued, must be positive, 0 to 9
--yearday, -t:      Year and day of year for plasmapause location (YYYYDDD)
--msec, -u:         Milliseconds into day, for plasmapause location
--geom_factor, -v:  Include the dipole geometric focusing factor? 1 for yes, 0 for no
----------------------------------------------------------- 
```
When using the modern mode (--mode=1), the file ```crres_clean.mat``` must be present in the current working directory.
The output is a text file, which contains two column vectors: time, and relative attenuation.

## dump_psd_models
This module writes the phase-space density model used in the modern implementation to a text file. Run it from the 
command line with the following arguments:

```

-------- Stanford Ray Tracer PSD Model Dumper ---------- 
  A program to export phase-space density models used
             by the Landau damping code.
 --out_file, -o:     Output file to use (default: dump.dat)
 --mode, -m:         Mode -- 0 for the legacy Suprathermal model,
                     1 for the CRRES model, and 2 for the hybrid model
 --AE, -a:           AE index - integer valued, 1, 2, 3.
 --Kp, -k:           Kp: Real valued, must be positive, 0 to 9
 --yearday, -t:      Year and day of year for plasmapause location (YYYYDDD)
 --msec, -u:         Milliseconds into day, for plasmapause location
 --MLT, -m:          MLT (Magnetic Local Time)
 --alpha, -p:        Electron pitch angle, in degrees
 --Lmin, -l:         Minimum L-shell to dump
 --Lmax, -L:         Maximum L-shell to dump
 --Lstep, -S:        L-shell step size
----------------------------------------------------------- 
```

The resulting output file is a space-separated list of entries, with each row containing: L-shell, energy, and f.
