#include <damping.h>
#include <gauss_legendre.h>

#include <complex>
#include <cmath>

using namespace std;
using namespace Eigen;

// extern "C" void sm_to_mag_d_(int* itime, double* x_in, double* x_out);
// extern "C" void cart_to_pol_d_(double* x_in, double* lat, double* lon, double* radial);

class psd_model; // Suprathermal electron distribution
class integrand; // The function to be integrated

// A port of Forrest's 3d damping code.
void damping_foust(rayF &ray, double Kp, double AE_level, int itime_in[2]) {
    double lat, lon, r;
    double L_sh;
    double mlt;

    Vector3d pos;
    Vector3d B0;
    Vector3d n_vec;
    Vector3d k;
    Vector3d Ns;
    Vector3d vgrel;
    Vector3d Bhat;
    Vector3d pos_prev;
    double pos_tmp[3]; 

    double kpar, kperp;
    double theta;
    double sin_th, cos_th, sin_th_sq, cos_th_sq;

    double n;

    double k_im, damp;
    double init_pwr = 1.0;   // Initial ray power

    double fs;
    double wce_h;
    double kmag;

    double R, L, P, S, D; // Stix parameters
    double a, b;          
    double B0mag;
    double wps2, whs;

    // results of integration
    double Di, ki;
    double ki_along_vg;
    // Integrand object
    integrand integ;

    // Resonance modes to integrate over
    // (0 = Landau damping, +-1 = cyclotron resonance)
    int m_low = -1;
    int m_hi  = 1;

    // Change this to an input, you doof
    // double AE_level = 3;

    // itime_in[0] = 2012045;          // yyyyddd
    // itime_in[1] = 1*(60*60*1000);   // time of day in msec     

    // n_steps = 500.0;
    // v_step = C/n_steps; //<- change back to this!

    // // Get ray launch location in mag dipole coordinates:
    // xin[0] = ray.pos[0][0];
    // xin[1] = ray.pos[0][1];
    // xin[2] = ray.pos[0][2];
    
    // // // Map to magnetic dipole coordinates
    // sm_to_mag_d_(itime_in, xin, xout);
    // cart_to_pol_d_(xout, &lat_init, &lon_init, &r_init);
    // // printf("pos size: %d\n",ray.pos[0].size());
    // r_init/= R_E;

    // printf("mag coords: %0.3f, %0.3f, %0.3f\n", R2D*lat_init, R2D*lon_init, r_init);


    // I don't understand what this is yet.
    // AN = AN_CM_DIST * pow( 10.0 , (12.0-(4.0*Q_DIST)) );

    // Initialize ray power with zeros
    // VectorXd ray_pwr = VectorXd::Zero(ray.time.size());

    // Initial power
    ray.damping = vector<double> (ray.time.size(), 0.0);
    ray.damping[0] = 1.0;

    // get L-shell of plasmapause, based on kp:
    double L_pp = kp_to_pp(Kp);

    cout << "L_plasmapause: " << L_pp << "\n";


    // ----------- Here's how we get the phase-space density model: ----------
    //  (equivalent to:   [n_fit, An_fit] = get_fit_params(L, MLT, AE_level, false);
    //                     fe = @(vperp, vpar) crres_polar_hybrid_psd(vperp, vpar, n_fit, An_fit, L, L_pp);    

    // (sorry for the hard-coded path here)
    char crres_data_file[100] = "/shared/users/asousa/software/damping/data/crres_clean.mat";
    
    psd_model psd;
    psd.initialize(crres_data_file);

    // Step through the ray:
    for (int ii=1; ii < ray.time.size(); ii++) {

        B0    = Map<VectorXd>(ray.B0[ii].data(), 3,1);
        pos   = Map<VectorXd>(ray.pos[ii].data(),3,1);
        pos_prev = Map<VectorXd>(ray.pos[ii-1].data(),3,1);
        n_vec = Map<VectorXd>(ray.n[ii].data(),3,1);
        Ns    = Map<VectorXd>(ray.Ns[ii].data(),3,1);
        vgrel = Map<VectorXd>(ray.vgrel[ii].data(),3,1);


        // Get local L-shell:
        lat = atan(pos[2]/sqrt(pow(pos[0],2) + pow(pos[1],2)));
        r = pos.norm();
        L_sh = r/pow(cos(lat),2)/R_E;


        // Get lat, lon, alt in geomagnetic:
        sm_to__d_(itime_in, pos.data(), pos_tmp);
        cart_to_pol_d_(pos_tmp, &lat, &lon, &r);

        // cout << "orig: " << r << ",\t" << lat << endl;
        // cout << "xfrm: " << r << ",\t" << lat << " lon: " << lon*R2D <<  endl;

        // Get MLT:
        mlt = fmod((atan2(pos[1], pos[0]) + PI)/(2*PI)*24, 24); //MLT in hours; 0 MLT is in -x direction
        // cout << "MLT orig: " << mlt << endl;

        mlt = MLT(itime_in, lon);
        // cout << "MLT (mine): " << mlt << endl;
        // cout << endl;

        // Set the current location parameters for the density model:
        L_pp = bulge(Kp, mlt); // Get plasmapause at this MLT
        // cout << "i: " << ii << "\tMLT: " << mlt<< "\tLpp: " << L_pp << endl;
        psd.set_params(L_sh, L_pp, mlt, AE_level);

        wce_h = Q_EL*B0.norm()/M_EL;

        k = n_vec*ray.w/C;
        kmag = k.norm();
        Bhat = B0.array()/B0.norm();
        kpar = k.dot(Bhat); //k.array()*Bhat.array();
        kperp = (k - kpar*Bhat).norm();

        B0mag = B0.norm();

        // ------- spatialdamping.m -----------
        // Theta is the angle between parallel and perpendicular K
        theta = atan2(kperp, kpar);

        // Some trig.
        sin_th = sin(theta);
        cos_th = cos(theta);
        sin_th_sq = pow(sin_th,2);
        cos_th_sq = pow(cos_th,2);


        n = n_vec.norm();
        // cout << "n (out): " << n << "\n";

        // ---------- Evaluate Stix parameters: ------
        wps2 = 0;
        R = 1.;
        L = 1.;
        P = 1.;

        for (int jj=0; jj < ray.Ns[ii].size(); jj++) {
            
            // Ns*(Qs^2)/(ms*Eps_0)
            wps2 = ray.Ns[ii][jj]*pow(ray.qs[jj],2) \
                   /(ray.ms[jj]*EPS0);
            // qB/m
            whs  = ray.qs[jj]*B0mag/ray.ms[jj];

            // Complex modification to whs -- for now, ignore. (8.19.2016)
            // wcs  = whs * ray.w/(ray.w + 1i*ray.nus[ii][jj]);

            R-= wps2/(ray.w*(ray.w + whs));
            L-= wps2/(ray.w*(ray.w - whs));
            P-= wps2/(ray.w*ray.w);
        }
        S = (R + L)/2.;
        D = (R - L)/2.;

        a = S*sin_th_sq + P*cos_th_sq;
        b = R*L*sin_th_sq + P*S*(1+cos_th_sq);


        // if (ii==1) {
        //             cout << "wps2: " << wps2 << " whs: " << whs << "\n";
        //             cout << "Stix Params [1]: ";
        //             cout << R << " ";
        //             cout << L << " ";
        //             cout << P << " ";
        //             cout << S << " ";
        //             cout << D << " ";
        //             cout << a << " ";
        //             cout << b << "\n";
        // }
        // --------------------------------------------

        // ---------- hot_dispersion_imag.m ------

        integ.initialize(psd, kperp, kpar, 
                        ray.w, n, m_low, m_hi, wce_h, 
                        R, L, P, S);

        // // Integrate it!
        Di = gauss_legendre(50, integrand_wrapper, (void*) &integ, 0., 1.);
        ki = -(ray.w/C)*(0.5)*(1./(4.*n*(2.*a*n*n-b))) * Di;
        // (ki is the output of spatialdamping.m)

        // ki_along_vg = ki*(k*vgrel(ii,:)')/(norm(k)*norm(vgrel(ii,:)));
        ki_along_vg = ki*(k.dot(vgrel)/(k.norm()*vgrel.norm()));

        double dist = (pos - pos_prev).norm();

        // cout << "dist: " << dist << "\n";
        // Power = previous power *(e^(-(distance)(damping)(2))).
        ray.damping[ii] = ray.damping[ii-1]*exp(-dist*ki_along_vg*2.0);        



    } // Step through ray


} // damping_foust


// Here's a working place to write an integrand. But how do we point it to
// a method of an object? Hmm.
double integrand_wrapper(double x, void* data) {

    integrand* integ = (integrand*)data;
    return integ->evaluate_t(x);

}
