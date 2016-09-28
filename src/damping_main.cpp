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
{
    map <int, rayF> raylist;
    // map <int, VectorXd> damplist;


    double x_in[3];
    double x_out[3];
    int itime;

    // char *fileName;
    string inpFileName;
    string outFileName;
    int mode;
    char damping_fileName[100];
    // ostringstream damping_fileName;
    FILE * outputFile;

    double AE_level;
    double Kp;

    // Default parameters:
    inpFileName = "four_adjacent.ray";
    outFileName = "damping.txt";
    mode = 1;     // 1 for foust, 0 for ngo
    AE_level = 3;
    Kp = 4;

    // Parse input arguments:
    int opt = 0;
    while ((opt = getopt(argc, argv, "i:o:m:a:k:")) != -1) {
        switch(opt) {
            case 'i':
            // input filename:
            // cout << "inp filename!\n";
                inpFileName = (string) optarg;
                break;
            case 'o':
            // output filename:
                outFileName = (string) optarg;
                break;
            case 'm':
            // damping mode:
                mode = atoi(optarg);
                break;
            case 'a':
                AE_level = strtod(optarg, NULL);
                break;
            case 'k':
                Kp = strtod(optarg, NULL);
                break;
            case '?':
                 printf("\nUnknown option: %s\n",opt);
            break;
        }
    }

    // Display current inputs:
    cout << "---- Input Parameters ----\n";
    cout << "input file: " << inpFileName << "\n";
    cout << "output file: " << outFileName << "\n";
    cout << "damping mode: " << mode << "\n";
    cout << "AE_level: " << AE_level << "\n";
    cout << "Kp: " << Kp << "\n";

    cout << "\n---- DAMPING ----\n";
    // Load the rayfile:

    raylist = read_rayfile(inpFileName);
    
    for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){
        printf("damping ray # %d\n",iter->first);

        // Calculate damping for the ray:
        // results are stored in ray.damping
        switch(mode) {
            case 0:
                damping_ngo(iter->second);   break;
            case 1:
                damping_foust(iter->second, Kp, AE_level); break;
        }
    }   


    // Print to file:
    write_rayfile(outFileName, raylist);

    // // Print to a file:
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter){

    //     sprintf(damping_fileName, "damping_%2.0f.txt",iter->second.w/(2*PI));
    //     outputFile = fopen(damping_fileName, "w");

    //     if (outputFile != NULL) {
    //         printf("Hey! Opened %s!\n",damping_fileName);

    //         for (int i=0; i < iter->second.damping.size(); i++) {
    //             fprintf(outputFile,"%0.4f\t%0.4f\n",iter->second.time[i], iter->second.damping[i]);
    //         }
    //         fclose(outputFile);

    //     } else {
    //         printf("Ugh, didn't open file!\n");
    //     }
    // }   



    // printf("x_in is: %g, %g, %g\n",x_in[0], x_in[1], x_in[2]);
    // sm_to_geo_d_(&itime, x_in, x_out);
    // printf("x_out is: %g, %g, %g\n",x_out[0], x_out[1], x_out[2]);
    // cart_to_pol_d_(&x_out, x_out[0], x_out[1], x_out[2]);
    // printf("x_out is: %g, %g, %g\n",x_out[0]*180./3.14, x_out[1]*180./3.14, x_out[2]/1000.);




    // Let's check some coordinate transforms!

    // double lat_in, lon_in, rad_in;
    // double tmp_in[3], tmp_out[3];
    // double tmp_out_2[3];

    // int itime_in[2];
    // int iyr, iday;
    // double ut, jd;


    // itime_in[0] = 2012001;          // yyyyddd
    // itime_in[1] = 1*(60*60*1000);   // time of day in msec

    // lat_in = 45*D2R;
    // lon_in = -180*D2R;
    // rad_in = 6371e3;

    // double lat_out, lon_out, rad_out;


    // pol_to_cart_d_(&lat_in, &lon_in, &rad_in, tmp_in);
    // geo_to_sm_d_(itime_in, tmp_in, tmp_out);
    // sm_to_geo_d_(itime_in, tmp_out, tmp_out_2);
    // cart_to_pol_d_(tmp_out_2, &lat_out, &lon_out, &rad_out);
    
    // lat_out = R2D*lat_out;
    // lon_out = R2D*lon_out;

    // printf("iyr: %d, iday: %d\n",iyr, iday);
    // printf("jd: %g\n",jd);
    // printf("tmp_out is: %g, %g, %g\n",tmp_out[0], tmp_out[1], tmp_out[2]);
    // printf("tmp_out_2 is: %g, %g, %g\n",tmp_out_2[0], tmp_out_2[1], tmp_out_2[2]);

    // printf("returned lat: %g, lon: %g, rad: %g\n",lat_out, lon_out, rad_out);



    // // Confirm we loaded everything
    // for(map<int,rayF>::iterator iter = raylist.begin(); iter != raylist.end(); ++iter)
    //     {
    //         float key = iter->first;
    //         rayF val = iter->second;
    //         printf("Ray number %g: freq: %g nspec: %g\n",key, val.w, val.nspec);

    //         // Print out all elements in a vector
    //         vector<float> vec = val.time;
    //         for (vector<float>::iterator it = vec.begin(); it != vec.end(); ++it)
    //             printf("%g ",*it);
    //         cout << "\n";

    //     }


    return 0; // Return statement.
} // Closing Main.

