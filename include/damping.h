// wipp.h
#ifndef damping_H
#define damping_H
#include <coord_transforms.h>

#include <Eigen/Core>
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
#include <getopt.h>

#include <integrand.h>
#include <psd_model.h>
#include <mat.h>
#include <consts.h>

using namespace std;

// Structure for holding an entire rayfile (yuge)
typedef struct rayfile_struct {
    int size;
    double w;                          // frequency (angular)
    double stopcond;                   // Stop condition
    double nspec;                      // Number of species in plasmasphere model
    // double ray_num;                 // Ray index number
    vector <double> time;            // Group time

    int iyr;
    int idoy;
    int isec;

    vector <vector <double> > pos;
    vector <vector <double> > vprel;
    vector <vector <double> > vgrel;
    vector <vector <double> > n;
    vector <vector <double> > B0;

    // Variable-length stuff (depending on number of constituents in model)
    vector <double> qs;    // species charge
    vector <double> ms;    // species mass
    vector <vector <double> > Ns;    // number density of species (m^-3)
    vector <vector <double> > nus;   // collision frequencies
    vector <double> damping;

} rayF;

// wipp_fileutils:
map<int, rayF> read_rayfile(string fileName);
void write_rayfile(string fileName, map <int, rayF> raylist);   // Write complete rayfile + damping
void write_damping(string fileName, map <int, rayF> raylist);   // Write damping standalone

// Damping (Ngo version):
// void damping_ngo(rayF &rayfile);
void damping_ngo(int itime_in[2], rayF &ray, bool include_geom_factor);


// Damping (Foust version):
void damping_foust(rayF &rayfile, double Kp, double AE_level, int itime_in[2], bool include_geom_factor);
double integrand_wrapper(double x, void* data);
double kp_to_pp(double kp);


// Math functions
double l2_norm(vector<double> u);
vector<double> scalar_multiply(vector<double> u, double v);
double dot_product(vector<double>u, vector<double>v);
vector<double> add(vector<double>u, vector<double> v);

// Helpers
void print_vector(vector<double> u);
void polyfit(const vector<double> &xv, const vector<double> &yv, vector<double> &coeff, int order);

// External Fortran function to calculate plasmapause vs Kp.
// Lifted from GCPM plasma model, GCPMv2.4
// extern "C" double pp_profile_d_(double* al, double* amlt, double* akp, double* a8);

double bulge(double kp, double mlt);

#endif
