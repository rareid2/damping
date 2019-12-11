#include <damping.h>

#include <sys/types.h>
#include <dirent.h>

#include "mat.h"    // Matlab file-access libraries


using namespace std;
using namespace Eigen;

void psd_model::initialize(char* inp_filename) {
    MATFile *pmat;
    DIR* dp;
    
    this->filename = inp_filename;
    const mwSize *dims;
    
        // Load CRRES data (cleaned up pls!)
        pmat = matOpen(inp_filename, "r");
        if (pmat == NULL) {
            cout << "Error opening file at " << inp_filename << "\n";
        } else {
            crres_data = matGetVariable(pmat, "crres");
            num_entries = mxGetNumberOfElements(crres_data);

            cout << "Loaded " << inp_filename << "\n";

            // dims = mxGetDimensions(crres_data);

            // cout << "ndims: " << mxGetNumberOfFields(crres_data) << "\n";
            // cout << dims[0] << " " << dims[1] << "\n";
        }
}

void psd_model::set_params(double L_sh_in, double L_pp_in, double MLT_in, double AE_level_in) {

    L_sh = L_sh_in;
    L_pp = L_pp_in;
    MLT  = MLT_in;
    AE_level = AE_level_in;

    vector<double> fit_params = CRRES_fit_params(L_sh, MLT, AE_level);
    n_fit = fit_params[0];
    An_fit= fit_params[1];
    // cout << "Fit parameters: L: " << L_sh << " MLT: " << MLT << " n: " << n_fit << ", An: " << An_fit << "\n";
    
}
void psd_model::replace_NaNs(mxArray* arr) {
    // Allegedly the CRRES data uses 1e12 to denote a NaN.
    /* Declare variables */ 
    size_t elements;
    mwSize j,cmplx;
    mwSize number_of_dims;
    mwSize nnz=0, count=0; 
    double *pr, *pi, *pind;
    mxComplexDouble *pc;    // Complex double -- upgrading to Matlab API 2019 
    const mwSize *dim_array;

    /* Get the data */
    // pc = mxGetComplexDoubles(arr);
    pr = (double *)mxGetDoubles(arr);
    // pr=(double *)mxGetPr(arr);
    // pi=(double *)mxGetPi(arr);
    // cmplx = ((pi==NULL) ? 0 : 1);

    number_of_dims=mxGetNumberOfDimensions(arr);
    elements=mxGetNumberOfElements(arr);
    // cout << " n_dims: " << number_of_dims;
    // cout << " elements: " << elements << "\n";

    for(int j=0;j<elements;j++) {
        if (pr[j]==1e12) {
            pr[j] = NAN;
            // cout << " nan @ " << j;
        }
    }
}

void psd_model::count_NaNs(mxArray* arr) {
    // Allegedly the CRRES data uses 1e12 to denote a NaN.
    /* Declare variables */ 
    size_t elements;
    mwSize j,cmplx;
    mwSize number_of_dims;
    mwSize nnz=0, count=0; 
    double *pr, *pi, *pind;
    const mwSize *dim_array;

    /* Get the data */
    // pr=(double *)mxGetPr(arr);
    // pi=(double *)mxGetPi(arr);
    // cmplx = ((pi==NULL) ? 0 : 1);
    pr = (double *)mxGetDoubles(arr);

    // cout << "complex? " << cmplx;

    number_of_dims=mxGetNumberOfDimensions(arr);
    elements=mxGetNumberOfElements(arr);
    // cout << " n_dims: " << number_of_dims;
    // cout << " elements: " << elements << "\n";
    int nanTotal = 0;
    for(int j=0;j<elements;j++) {
        if (mxIsNaN(pr[j])) {
            nanTotal++;
        }
    }
    cout << "Total NaNs: " << nanTotal << "\n";
}



vector<double> psd_model::CRRES_fit_params(double L, double MLT, double AE_level) {
    // Loop through each of the CRRES entries, and find
    // the closest values for E and J_perp.
    // 
    // This is a port of "get_fit_params.m" by dgolden.

    int row, col;
    mxArray *L_mx;
    mxArray *MLT_mx;
    mxArray *Jp_mx;
    mxArray *E_mx;
    mxArray *AE_mx;


    double *MLTd;
    double *Ld;
    double *Jpd;
    MatrixXd MLT_data, L_data, Jp_data;
    double *E_d;
    double *AE_d;

    const mwSize *dims;

    // VectorXd log_Jp(num_entries);
    vector<double> log_Jp;
    vector<double> log_E;

    vector<double> fit_params;
    vector<double> p;   // polyfit coefficients
    double J0;
    double m;

    double AE_target = round(AE_level);
    // cout << "AE_target: " << AE_target << "\n";

    // AE is either 1, 2, or 3 (after rounding)
    if ( (AE_target < 1) || (AE_target > 3) ) {
        cout << "AE_target out of range!\n";
        return fit_params;
    }

    for (int kk = 0; kk < num_entries; kk++) {
        // Get mxArrays from the structure
        MLT_mx = mxGetField(crres_data, kk, "MLT");         // MLT <mxn>
        L_mx   = mxGetField(crres_data, kk, "L");           // L-shell <mxn>
        Jp_mx  = mxGetField(crres_data, kk, "J_int");       // Jperp, interpolated <mxn>
        E_mx   = mxGetField(crres_data, kk, "E");           // Energy (keV) <double>
        AE_mx  = mxGetField(crres_data, kk, "AE");          // AE level     <double>

        // Get pointers to the sweet data within
        // MLTd = mxGetPr(MLT_mx);
        // Ld   = mxGetPr(L_mx);
        // Jpd  = mxGetPr(Jp_mx);
        // E_d  = mxGetPr(E_mx);
        // AE_d = mxGetPr(AE_mx);

        MLTd = mxGetDoubles(MLT_mx);
        Ld   = mxGetDoubles(L_mx);
        Jpd  = mxGetDoubles(Jp_mx);
        E_d  = mxGetDoubles(E_mx);
        AE_d = mxGetDoubles(AE_mx);

        if (AE_d[0] == AE_target) {

            // cout << "E: " << E_d[0] << " AE: "<< AE_d[0] << "\n";
            // get dimensions:
            dims = mxGetDimensions(MLT_mx);

            // Map to Eigen matrices
            MLT_data = Map<MatrixXd>(MLTd,dims[0], dims[1]);
            L_data   = Map<MatrixXd>(Ld,  dims[0], dims[1]);
            Jp_data  = Map<MatrixXd>(Jpd, dims[0], dims[1]);

            int number_of_dims=mxGetNumberOfDimensions(MLT_mx);
            int elements=mxGetNumberOfElements(MLT_mx);

            // Find the closest value to MLT in the first row of MLT_data:
            // (row, col indexes are probably constant for each run thru the CRRES data,
            // but let's leave this here for now unless it's obnoxiously slow)
            (MLT_data.row(0).array() - MLT).abs().minCoeff(&row);    

            // Same search for L vector:
            (L_data.col(0).array() - L).abs().minCoeff(&col);
           

            log_Jp.push_back(log(Jp_data(col, row)));
            log_E.push_back(log(E_d[0]));

        }   // AE level
    }   // Loop over each entry in crres_data
    
    // cout << "Log_E: \n";
    // print_vector(log_E);
    // cout << "Log_Jperp: \n";
    // print_vector(log_Jp);


    // Linear fit:
    polyfit(log_E, log_Jp, p, 1);
    
    // cout << "Polyfit coeffs: \n";
    // print_vector(p);

    J0 = exp(p[0]);
    m  = -1.0*p[1];

    // cout << "J0: " << J0 << " m: " << m << "\n";
    n = 2*m + 2;
    // Watch out here! M_EL is our electron mass in consts.h,
    // but there's an "M_E" defined somewhere too.
    An = 2*J0/pow((0.5*(6.25e11)*M_EL), m-1);
    
    fit_params.push_back(n);
    fit_params.push_back(An);

    // cout << "Fit Params: \n";
    // print_vector(fit_params);

    return fit_params;
}

double psd_model::suprathermal(double vperp, double vpar) {
    // Suprathermal distribution.
    // This is what we'll use inside the plasmasphere.

    double v, f;

    const double a = 4.9e5;
    const double b = 8.3e14;
    const double c = 5.4e23;
    // % Just a crutch to avoid the singularity
    const double v0=1;

    // % Convert to cm/s
    v = 100.*sqrt(vperp*vperp + vpar*vpar + v0);
    f = (a/pow(v,4) -b/pow(v,5) + c/pow(v,6));

    // % Convert to s^3/m^6 from s^3/cm^6 
    f = f*pow(100.,6);
    return f;
}

// double psd_model::crres_psd(double vperp, double vpar, double n, double An) {
double psd_model::crres_psd(double vperp, double vpar) {
    // Suprathermal distribution, model-driven.
    // This is what we'll use outside the plasmasphere.

    // % Get phase space density using CRRES suprathermal fluxes from Bortnik
    // % [2007]
    // % 
    // % vperp and vpar in m/s

    // % By Daniel Golden (dgolden1 at stanford dot edu) May 2010
    // % $Id: crres_psd.m 1205 2011-03-04 23:15:40Z dgolden $

    double v, f;

    // % Just a crutch to avoid the singularity
    const double v0=1;

    // % Convert to cm/s
    v = 100.*sqrt(vperp*vperp + vpar*vpar + v0);

    f = (this->An_fit)/pow(v, this->n_fit);

    // % Convert to s^3/m^6 from s^3/cm^6 
    f = f*pow(100.,6);
    return f;
}

// double psd_model::hybrid_psd(double vperp, double vpar, double n_fit, double An_fit, double L, double L_pp) {
double psd_model::hybrid_psd(double vperp, double vpar) {
    double f;
    double f_polar, f_crres;
    double w_polar, w_crres;

    if (L_pp - L_sh > 1) {
        // cout << "inside:\n";
        // Way inside plasmasphere:
        f = this->suprathermal(vperp, vpar);
    } else if (L_sh - L_pp > 1) {
        // cout << "outside:\n";
        // Way outside plasmasphere:
        // f = this->crres_psd(vperp, vpar, n_fit, An_fit);
        f = this->crres_psd(vperp, vpar);

    } else {
        // cout << "hybrid:\n";
        // Blend between the two:
        f_polar = this->suprathermal(vperp, vpar);
        // f_crres = this->crres_psd(vperp, vpar, n_fit, An_fit);
        f_crres = this->crres_psd(vperp, vpar);

        // cout << "f_polar: " << f_polar << " f_crres: " << f_crres << "\n";
        w_polar = exp(5*(L_pp - L_sh))/(1 + exp(5*(L_pp - L_sh))); // higher weight inside plasmasphere
        w_crres = exp(5*(L_sh - L_pp))/(1 + exp(5*(L_sh - L_pp))); // higher weight outside plasmasphere

        // Find centroid in log space
        f = exp( (log(f_polar)*w_polar + log(f_crres)*w_crres)
                         /(w_polar + w_crres) );
    }

    return f;
}


