#ifndef Utils
#define Utils

// This is from Byukling Kayanti Formula (6.3)
double Lambda(double x, double y, double z)
{
        return (x - y - z) * (x - y - z) - 4 * y * z;
}

// From Byukling Kayanti Formula (5.14) Page 86
double T_min(double ma_2, double mb_2, double m1_2, double m2_2, double s) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
        return ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) - sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2)));
}

double T_max(double ma_2, double mb_2, double m1_2, double m2_2, double s)
{
        return ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) + sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2)));
}

double T_min_function(double *x, double *par) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
        const double Mp = 0.9383;    // GeV
        const double M_JPsi = 3.097; // GeV
        double Eg = x[0];
        double ma_2 = 0.;
        double mb_2 = Mp * Mp;
        double m1_2 = M_JPsi * M_JPsi;
        double m2_2 = Mp * Mp;
        double s = Mp * Mp + 2 * Mp * Eg;
        return -1 * (ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) - sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))));
}

double T_max_function(double *x, double *par)
{
        const double Mp = 0.9383;    // GeV
        const double M_JPsi = 3.097; // GeV
        double Eg = x[0];
        double ma_2 = 0.;
        double mb_2 = Mp * Mp;
        double m1_2 = M_JPsi * M_JPsi;
        double m2_2 = Mp * Mp;
        double s = Mp * Mp + 2 * Mp * Eg;
        return -1 * (ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) + sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))));
}

double Bin_limit_T_min_function(double *x, double *par) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
        const double Mp = 0.9383;    // GeV
        const double M_JPsi = 3.097; // GeV
        double Eg = x[0];
        double ma_2 = 0.;
        double mb_2 = Mp * Mp;
        double m1_2 = M_JPsi * M_JPsi;
        double m2_2 = Mp * Mp;
        double s = Mp * Mp + 2 * Mp * Eg;
        double t_min = -1 * (ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) - sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))));
        return TMath::Max(t_min, par[0]);
}

double Bin_limit_T_max_function(double *x, double *par)
{
        const double Mp = 0.9383;    // GeV
        const double M_JPsi = 3.097; // GeV
        double Eg = x[0];
        double ma_2 = 0.;
        double mb_2 = Mp * Mp;
        double m1_2 = M_JPsi * M_JPsi;
        double m2_2 = Mp * Mp;
        double s = Mp * Mp + 2 * Mp * Eg;
        double t_max = -1 * (ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) + sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))));
        return TMath::Min(t_max, par[0]);
}

double bin_volume_correction(double min_bin_t, double max_bin_t, double min_bin_Eg, double max_bin_Eg, bool debug, TString name)
{
        TF1 *f_Bin_limit_T_min_function = new TF1("f_Bin_limit_T_min_function", Bin_limit_T_min_function, min_bin_Eg, max_bin_Eg, 1);
        TF1 *f_Bin_limit_T_max_function = new TF1("f_Bin_limit_T_max_function", Bin_limit_T_max_function, min_bin_Eg, max_bin_Eg, 1);

        // TH2F *histogram_integration = new TH2F("histogram_integration", "histogram_integration", min_bin_Eg, max_bin_Eg, 0.0, 8.0);

        f_Bin_limit_T_min_function->SetParameters(min_bin_t, 0.0);
        f_Bin_limit_T_max_function->SetParameters(max_bin_t, 0.0);

        TRandom3 randGen;

        // Number of random points for integration
        int numPoints = 100000;

        // Calculate the area between the two functions using Monte Carlo integration
        double sum = 0.0;
        for (int i = 0; i < numPoints; ++i)
        {
                double x = randGen.Uniform(min_bin_Eg, max_bin_Eg);
                double y_min = f_Bin_limit_T_min_function->Eval(x);
                double y_max = f_Bin_limit_T_max_function->Eval(x);
                sum += (y_max > y_min) ? std::abs(y_max - y_min) : 0.0;
        }

        // Scale the sum to get the estimated area
        double area = sum * (max_bin_Eg - min_bin_Eg) / numPoints;
        cout << " " << endl;
        cout << "area 1 " << area << endl;
        cout << "area 2 " << ((max_bin_t - min_bin_t) * (max_bin_Eg - min_bin_Eg)) << endl;

        if (debug)
        {
                TCanvas *canvas = new TCanvas("canvas", "Two Function Plot", 800, 600);
                canvas->cd();
                f_Bin_limit_T_min_function->SetLineColor(kBlue);
                f_Bin_limit_T_min_function->SetMinimum(min_bin_t);
                f_Bin_limit_T_min_function->SetMaximum(max_bin_t);
                f_Bin_limit_T_max_function->SetLineColor(kRed);
                f_Bin_limit_T_min_function->Draw();
                f_Bin_limit_T_max_function->Draw("SAME");
                canvas->SaveAs(name + ".pdf");
        }

        return area / ((max_bin_t - min_bin_t) * (max_bin_Eg - min_bin_Eg));
}

double JPsi_sigm_int_2g(double *xx, double *par)
{
        const double Mp = 0.9383;    // GeV
        const double M_JPsi = 3.097; // GeV
        const double PI = 3.14159;
        const double Fermi = 1.e-13; // cms

        const double R = 1. * Fermi;
        double N_2g = 1.;
        // const double t_slope = 1.13;

        double Eg = xx[0];
        N_2g = par[1];
        double t_slope = par[2];
        // double F_2g = TMath::Exp(t_slope*tM);
        double s = Mp * Mp + 2 * Mp * Eg;
        double x = (2 * Mp * M_JPsi + M_JPsi * M_JPsi) / (s - Mp * Mp);
        double Q2 = M_JPsi * M_JPsi;

        double nue = 1 / (16 * PI * (s - Mp * Mp) * (s - Mp * Mp));
        double t_min = T_min(0., Mp * Mp, Q2, Mp * Mp, s);
        double t_max = T_max(0., Mp * Mp, Q2, Mp * Mp, s);

        double Sigm_int = N_2g * nue * TMath::Power(1 - x, 2) / (R * R * M_JPsi * M_JPsi) * (1. / t_slope) * (TMath::Exp(t_slope * t_min) - TMath::Exp(t_slope * t_max)) * TMath::Power(s - Mp * Mp, 2);
        return Sigm_int;
}

double Fit_diff_CS_t(double *x, double *par)
{
        double t = -1.0*x[0];
        double formula_t = par[0]/TMath::Power((1.-(t)/TMath::Power(par[1],2)), 4);
        return formula_t;
}

//////////////////// Glue X and Hall C results ////////////////////////////////
/*float EgGluex[18] = {8.3, 8.5, 8.7, 8.8, 9.0, 9.2, 9.37, 9.55, 9.7, 9.9, 10.1, 10.3, 10.5, 10.6, 10.8, 11.0, 11.2, 11.35};
float xsGluex[18] = {0.044, 0.137, 0.252, 0.330, 0.226, 0.188, 0.510, 0.786, 0.518, 0.703, 0.795, 0.813, 1.376, 0.969, 1.142, 1.092, 1.563, 1.847};
float sxsGluex[18] = {0.028, 0.025, 0.028, 0.014, 0.056, 0.021, 0.022, 0.068, 0.021, 0.116, 0.073, 0.046, 0.066, 0.116, 0.051, 0.031, 0.130, 0.027};
float xeGluex[18] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};*/

float EgGluex[18] = {8.29, 8.47, 8.65, 8.83, 9.01, 9.19, 9.37, 9.55, 9.73, 9.91, 10.09, 10.27, 10.45, 10.63, 10.81, 10.99, 11.17, 11.35};
float xsGluex[18] = {0.043, 0.136, 0.249, 0.326, 0.206, 0.200, 0.489, 0.710, 0.507, 0.683, 0.829, 0.848, 1.321, 0.981, 1.151, 1.114, 1.594, 1.791};
float sxsGluex[18] = {0.012, 0.022, 0.029, 0.048, 0.059, 0.060, 0.087, 0.134, 0.08, 0.100, 0.119, 0.123, 0.193, 0.134, 0.140, 0.126, 0.208, 0.344};
float xeGluex[18] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};


float EgGluex_old[10] = {8.38, 8.74, 9.15, 9.46, 9.82, 10.18, 10.54, 10.9, 11.26, 11.62};
float xsGluex_old[10] = {0.116, 0.343, 0.313, 0.835, 0.868, 0.949, 1.383, 1.274, 2.158, 3.245};
float sxsGluex_old[10] = {0.031, 0.067, 0.127, 0.194, 0.196, 0.187, 0.284, 0.206, 0.421, 0.928};
float xeGluex_old[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

float Ei[10] = {9.175, 9.325, 9.475, 9.625, 9.775, 9.925, 10.075, 10.225, 10.375, 10.525};
float xstot[10] = {0.31291667, 0.47101143, 0.5038489, 0.59467919, 0.4996794, 0.55587942, 0.72679993, 0.72072183, 0.70057103, 0.79884448};
float sxstot[10] = {0.19454057, 0.14799775, 0.15983122, 0.16940087, 0.22375881, 0.23247598, 0.21068647, 0.2021271, 0.24060015, 0.28827128};
float xetot[18] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

float E_Richard[8] = {8.66983, 9.10284, 9.30166, 9.50429, 9.7028, 9.89812, 10.0908, 10.3569};
float xstot_Richard[8] = {0.022847, 0.208625, 0.166128, 0.235696, 0.219467, 0.215098, 0.375233, 0.399482};
float s_xstot_Richard[8] = {0.0110167, 0.0545851, 0.0522961, 0.0647822, 0.0689682, 0.0657605, 0.0894844, 0.182945};
float s_E_Richard[8] = {0.234984, 0.0559321, 0.0562257, 0.0578705, 0.0576846, 0.059219, 0.0584203, 0.100347};

float t_GluX_Bin2[9] = {0.69, 0.87, 1.21, 1.71, 2.24, 2.97, 3.89, 5.06, 6.37};
float xs_t_GluX_Bin2[9] = {0.813, 0.499, 0.401, 0.231, 0.120, 0.075, 0.026, 0.019, 0.009};
float e_xs_t_GluX_Bin2[9] = {0.088, 0.061, 0.037, 0.027, 0.021, 0.011, 0.008, 0.005, 0.004};
float e_t_up_GluX_Bin2[9] = {0.08, 0.13, 0.29, 0.29, 0.26, 0.53, 0.61, 0.69, 1.73};
float e_t_down_GluX_Bin2[9] = {0.2, 0.1, 0.21, 0.21, 0.24, 0.47, 0.39, 0.56, 0.62};
float e_t_up_GluX_Bin2_zeros[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
float e_t_down_GluX_Bin2_zeros[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

float t_GluX_Bin1[7] = {0.92, 1.25, 1.72, 2.24, 2.94, 3.92, 4.95};
float xs_t_GluX_Bin1[7] = {0.313, 0.170, 0.097, 0.045, 0.018, 0.030, 0.033};
float e_xs_t_GluX_Bin1[7] = {0.092, 0.018, 0.010, 0.007, 0.003, 0.006, 0.013};
float e_t_up_GluX_Bin1[7] = {0.08, 0.25, 0.28, 0.26, 0.56, 0.58, 0.80};
float e_t_down_GluX_Bin1[7] = {0.15, 0.25, 0.22, 0.24, 0.44, 0.42, 0.45};
float e_t_up_GluX_Bin1_zeros[7] = {0., 0., 0., 0., 0., 0., 0.};
float e_t_down_GluX_Bin1_zeros[7] = {0., 0., 0., 0., 0., 0., 0.};

//double w_c_from_BG_estimation = 0.69;
double w_c_from_BG_estimation = 0.77;
double convert_Gev_to_fm = 0.1973;
float Branching_ratio = 0.06;

#endif
