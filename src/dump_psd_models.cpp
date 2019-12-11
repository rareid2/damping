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
// Version 1.0  9.2016  Austin Sousa  (asousa@stanford.edu)

// Input parameters:
//     -i:     Input file to use  (default: "input.ray")
//     -o:     Output file to use (default: "output.ray")
//     -m:     Mode -- 1 for modern implementation, 0 for the legacy
//             damping code from the 2d raytracer.
//     -a:     AE index - integer valued, 1, 2, 3.
//     -k:     Kp: Real valued, must be positive.

{

    map <int, rayF> raylist;


    double x_in[3];
    double x_out[3];
    int itime_in[2];

    // char *fileName;
    string inpFileName;
    string outFileName;
    int mode;
    char damping_fileName[100];
    // ostringstream damping_fileName;
    FILE * outputFile;

    double AE_level;
    double Kp;
    double L_sh;
    double MLT;
    double L_min;
    double L_max;
    double L_step;
    double E, E_min, E_max, dE;
    double f;
    double v_tot;

    double vperp, vpar;
    double alpha;


    // Default parameters:
    
    outFileName = "dump.dat";
    mode = 1;     // 1 for foust, 0 for ngo
    AE_level = 3;
    Kp = 4;
    L_min = 1.1;
    L_max = 5.0;
    L_step= 0.05;

    E_min = 10;
    E_max = 1e8;
    dE = 0.05;
    MLT = 0;
    alpha = D2R* 45;


    // Parse input arguments:
    static struct option long_options[] = 
    {
        {"out_file",    required_argument,  0,  'o'},
        {"mode",        required_argument,  0,  'm'},
        {"AE",          required_argument,  0,  'a'},
        {"Kp",          required_argument,  0,  'k'},
        {"yearday",     required_argument,  0,  't'},
        {"msec",        required_argument,  0,  'u'},
        // {"L_sh",        required_argument,  0,  'L'},
        {"MLT",         required_argument,  0,  'M'},
        {"alpha",       required_argument,  0,  'p'},
        {"Lmin",        optional_argument,  0,  'l'},
        {"Lmax",        optional_argument,  0,  'L'},
        {"Lstep",       optional_argument,  0,  'S'},
        {"Emin",        optional_argument,  0,  'e'},
        {"Emax",        optional_argument,  0,  'E'},
        {"Estep",       optional_argument,  0,  's'},
        {0, 0, 0, 0}
    };


    int opt = 0;
    int opt_index = 0;

    while (( opt = getopt_long (argc, argv, "o:m:a:k:t:u:", long_options, &opt_index)) != -1) {
        switch(opt) {
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
            case 'l':   // Lmin
                L_min = strtod(optarg, NULL);               break;
            case 'L':   // Lmax
                L_max = strtod(optarg, NULL);               break;
            case 'S':   // Lstep
                L_step = strtod(optarg, NULL);              break;
            case 'e':   // Emin
                E_min = strtod(optarg, NULL);               break;
            case 'E':   // Emax
                E_max = strtod(optarg, NULL);               break;
            case 's':   // Estep
                dE = strtod(optarg, NULL);                  break;                
            case 'M':   // msec of day
                MLT = strtod(optarg, NULL);                 break;
            case 'p':   // msec of day
                alpha = D2R*strtod(optarg, NULL);           break;
            case '?':
                 // printf("\nUnknown option: %s\n",opt);
            break;
        }
    }




        // ---------- Write a welcome message + help string: ----------
    cout << "\n\n";
    cout << "-------- Stanford Ray Tracer PSD Model Dumper ---------- \n";
    cout << "  A program to export phase-space density models used\n";
    cout << "             by the Landau damping code.\n";
    cout <<" --out_file, -o:     Output file to use (default: ""dump.dat"")\n";
    cout <<" --mode, -m:         Mode -- 0 for the legacy Suprathermal model\n";
    cout <<"                     1 for the CRRES model, and 2 for the hybrid model\n";
    cout <<" --AE, -a:           AE index - integer valued, 1, 2, 3.\n";
    cout <<" --Kp, -k:           Kp: Real valued, must be positive, 0 to 9\n";
    cout <<" --yearday, -t:      Year and day of year for plasmapause location (YYYYDDD)\n";
    cout <<" --msec, -u:         Milliseconds into day, for plasmapause location\n";
    cout <<" --MLT, -m:          MLT (Magnetic Local Time)\n";
    cout <<" --alpha, -p:        Electron pitch angle, in degrees\n";
    cout <<" --Lmin, -l:         Minimum L-shell to dump\n";
    cout <<" --Lmax, -L:         Maximum L-shell to dump\n";
    cout <<" --Lstep, -S:        L-shell step size\n";

    cout << "----------------------------------------------------------- \n";

    cout <<"\n\n";

    // Display current inputs:
    cout << "---- Input Parameters ----\n";
    cout << "output file: " << outFileName << "\n";
    cout << "damping mode: " << mode << "\n";
    cout << "AE_level: " << AE_level << "\n";
    cout << "Kp: " << Kp << "\n";
    cout << "Alpha: " << R2D*alpha << "\n";
    cout << "L_min: " << L_min << " L_max: " << L_max << " L_step: " << L_step <<"\n";
    cout << "E_min: " << E_min << " E_max: " << E_max << " dE: " << dE << "\n";
    // cout << "E_min: " << E_MIN <<" eV E_max: " << E_MAX << " eV num E: " << NUM_E << "\n";

    cout << "\n---- PSD Model Dump ----\n";
    
    // (hard-coded path, sorry)
    char crres_data_file[100] = "crres_clean.mat";

    // get pitch angle in radians
    alpha = D2R*alpha;

    // Initialize PSD model
    psd_model psd;
    psd.initialize(crres_data_file);

    // get plasmapause location
    double L_pp = kp_to_pp(Kp);

    cout << "for Kp = "<< Kp <<", plasmapause is at L=" << L_pp << "\n";


    outputFile = fopen(outFileName.c_str(), "w");    

    for (L_sh = L_min; L_sh < L_max; L_sh+=L_step) {
        cout << "computing for L = " << L_sh << "\n";
        // for (E = E_EXP_BOT; E < E_EXP_TOP; E+=2.0*DE_EXP) {
        for (E = log10(E_min); E < log10(E_max); E+=dE) {

            psd.set_params(L_sh, L_pp, MLT, AE_level);
            v_tot = C*sqrt(1 - pow( (E_EL/(E_EL+pow(10, E))) ,2) );

            // pitch angle -- angle between perp and par velocities
            vperp = v_tot*sin(alpha);
            vpar  = v_tot*cos(alpha);
            switch(mode) {
            
            case 0:    // Bell 2002 model ('suprathermal')
                f = psd.suprathermal(vperp, vpar);
                break;
            case 1:
                f = psd.crres_psd(vperp, vpar);
                break;
            case 2:
                f = psd.hybrid_psd(vperp, vpar);
                break;
            default:
                cout << "oops";
                break;
            }
            // Write the computed value to the file:
            fprintf(outputFile, "%g %g %g\n",L_sh, E, f);

        // Uncomment this for verbose output   
        // cout << "L: " << L_sh << " E: " << E << " f: " << f << "\n";
        }
    }
    fclose(outputFile);
    return 0; // Return statement.
} // Closing Main.

