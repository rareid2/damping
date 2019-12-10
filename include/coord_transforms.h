// wipp.h
#ifndef coord_xforms_H
#define coord_xforms_H
#include <Eigen/Core>
#include <Eigen/Dense>  // Cross product lives here

#include <algorithm>    // std::next_permutation, std::sort
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <map>
// #include <ctime>
#include <getopt.h>
#include <consts.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ftw.h>

using namespace std;
using namespace Eigen;

// ---- Coordinate transforms ----
extern "C" void geo2mag1_(int* iyr, double* xGEO, double* xMAG);

// ----- libxformd (Coordinate transform library used in the raytracer) -----
extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);
extern "C" void geo_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_sm_d_(int* itime, double* x_in, double* x_out);

extern "C" void geo_to_mag_d_(int* itime, double* x_in, double* x_out);
extern "C" void mag_to_geo_d_(int* itime, double* x_in, double* x_out);

// (in radians)
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radius);
extern "C" void pol_to_cart_d_(double* lat, double* lon, double* radius, double* x_out);

extern "C" void sm_to_gsm_d_(int* itime, double* x_in, double* x_out);
extern "C" void gsm_to_sm_d_(int* itime, double* x_in, double* x_out);

// ---- My own transforms ----
// In-place cartesian / polar transforms. 
void carsph(double x[3]); 
void sphcar(double x[3]); 
void cardeg(double x[3]);
void degcar(double x[3]);

void carsph(double x[3], double x_out[3]); 
void sphcar(double x[3], double x_out[3]); 
void cardeg(double x[3], double x_out[3]); 
void degcar(double x[3], double x_out[3]); 

// In-place mapping of a data field between cartesian / polar frames.
void transform_data_sph2car(double lat, double lon, double d_in[3], double d_out[3]);
void transform_data_car2sph(double lat, double lon, double d_in[3], double d_out[3]);
void transform_data_geo2mag(int itime_in[2], double d_in[3], double d_out[3]);
void transform_data_mag2geo(int itime_in[2], double d_in[3], double d_out[3]);

double MLT(int itime[2], double lon);


#endif