#include <damping.h>

using namespace std;
using namespace Eigen;


// External functions we'll use (libxformd for coordinate transforms)
extern "C" void sm_to_geo_d_(int* itime, double* x_in, double* x_out);
extern "C" void geo_to_sm_d_(int* itime, double* x_in, double* x_out);
extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radius);
extern "C" void pol_to_cart_d_(double* lat, double* lon, double* radius, double* x_out);
// extern "C" double t0_d_(int* itime, int* iyr, int* iday, double* ut); // date to julian centuries
 
int main(int argc, char *argv[]) 
// Calculates wave power damping for the Stanford 3D raytracer.
// A port of Forrest Foust's Matlab-based implementation.
// Input parameters:
//     --inp_file, -i:     Input file to use  (default: "input.ray")
//     --out_file, -o:     Output file to use (default: "output.ray")
//     --mode, -m:         Mode -- 1 for modern implementation, 0 for the legacy
//                         damping code from the 2d raytracer.
//     --AE, -a:           AE index - integer valued, 1, 2, 3.
//     --Kp, -k:           Kp: Real valued, must be positive.
//     --yearday, -t:      Year and day of year for plasmapause location (YYYYDDD)
//     --msec, -u:         Milliseconds into day, for plasmapause location
//     --geom_factor, -v:  Include the dipole geometric focusing factor? 1 for yes, 0 for no

// Version 1.0  9.2016  Austin Sousa  (asousa@stanford.edu)
//        -- initial port
// Version 1.1  12.2019 Austin Sousa  (Austin.Sousa@colorado.edu)
//        -- Some code cleanup
//        -- updated to use Matlab 2018+ for the file loading
//        -- briefly confirmed that the output matches the original Matlab code
//

{
    map <int, rayF> raylist;

    double x_in[3];
    double x_out[3];
    int itime_in[2];

    string inpFileName;
    string outFileName;
    int mode;
    char damping_fileName[100];
    FILE * outputFile;

    double AE_level;
    double Kp;
    bool include_geom_factor;

    // ---------- Default parameters: ----------
    itime_in[0] = 2010001;
    itime_in[1] = 0;
    inpFileName = "input.ray";
    outFileName = "output.damp";
    mode = 1;     // 1 for foust, 0 for ngo
    AE_level = 3;
    Kp = 4;
    include_geom_factor = false;

    // ---------- Parse input arguments: ----------
    static struct option long_options[] = 
    {
        {"inp_file",    required_argument,  0,  'i'},
        {"out_file",    required_argument,  0,  'o'},
        {"mode",        required_argument,  0,  'm'},
        {"AE",          required_argument,  0,  'a'},
        {"Kp",          required_argument,  0,  'k'},
        {"yearday",     required_argument,  0,  't'},
        {"msec",        required_argument,  0,  'u'},
        {"geom_factor", required_argument,  0,  'v'},
        {0, 0, 0, 0}
    };

    int opt = 0;
    int opt_index = 0;

    while (( opt = getopt_long (argc, argv, "i:o:m:a:k:t:u:", long_options, &opt_index)) != -1) {
        switch(opt) {
            case 'i':   // input filename:
                inpFileName = (string) optarg;              break;
            case 'o':   // output filename:
                outFileName = (string) optarg;              break;
            case 'm':   // damping mode:
                mode = atoi(optarg);                        break;
            case 'a':   // AE
                AE_level = strtod(optarg, NULL);            break;
            case 'k':   // Kp
                Kp = strtod(optarg, NULL);                  break;
            case 't':   // yearday
                itime_in[0] = atoi(optarg);                 break;
            case 'u':   // msec of day
                itime_in[1] = atoi(optarg);                 break;
            case 'v':   // include geometric factor?
                include_geom_factor = (atoi(optarg)==1);    break;    
            case '?':
                 // printf("\nUnknown option: %s\n",opt);
            break;
        }
    }

    // ---------- Write a welcome message + help string: ----------
    cout << "\n\n";
    cout << "-------- Stanford Ray Tracer Landau Damping Code ---------- \n";
    cout <<"--inp_file, -i:     Input file to use  (default: ""input.ray"")\n";
    cout <<"--out_file, -o:     Output file to use (default: ""output.ray"")\n";
    cout <<"--mode, -m:         Mode -- 1 for modern implementation, 0 for the legacy\n";
    cout <<"                    damping code from the 2d raytracer.\n";
    cout <<"--AE, -a:           AE index - integer valued, 1, 2, 3.\n";
    cout <<"--Kp, -k:           Kp: Real valued, must be positive, 0 to 9\n";
    cout <<"--yearday, -t:      Year and day of year for plasmapause location (YYYYDDD)\n";
    cout <<"--msec, -u:         Milliseconds into day, for plasmapause location\n";
    cout <<"--geom_factor, -v:  Include the dipole geometric focusing factor? 1 for yes, 0 for no\n";
    cout << "----------------------------------------------------------- \n";

    cout <<"\n\n";

    // ---------- Display current inputs: ----------
    cout << "---- Input Parameters ----\n";
    cout << "input file: " << inpFileName << "\n";
    cout << "output file: " << outFileName << "\n";
    cout << "damping mode: " << mode << "\n";
    cout << "AE_level: " << AE_level << "\n";
    cout << "Kp: " << Kp << "\n";
    cout << "Geometric factor? " << include_geom_factor <<"\n";
    cout << "\n---- DAMPING ----\n";

    // ---------- Load the rayfile: ----------
    raylist = read_rayfile(inpFileName);
    
    // Iterate over each entry in the ray file:
    // (iter->second is the current value of the iterator)
    for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
        printf("damping ray # %d\n",iter->first);

        // ---------- Calculate damping for the ray: ----------
        switch(mode) {
            case 0:
                damping_ngo(itime_in, iter->second, include_geom_factor);   break;
            case 1:
                damping_foust(iter->second, Kp, AE_level, itime_in, include_geom_factor); break;
        }
    }   


    // ---------- Print to file: ----------
    // write_rayfile(outFileName, raylist);     // Annotate whole rayfile
    write_damping(outFileName, raylist);        // Write the damping file separately
    
    return 0; // Return statement.
} // Closing Main.

