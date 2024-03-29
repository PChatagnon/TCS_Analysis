#ifndef Sample_Class
#define Sample_Class

class Sample
{
public:
        TString variable;     // Name of the variable
       
        float cut_acc = 0.05;
        float cut_acc_error = 0.5;

        TString outputfolder = "Plots_Pass2"; // outputfolder for plots and graphs

        /*Sample(TString input_variable, int input_bin_number, float *input_binning_array)
        {
                variable = input_variable;
                binning_array = input_binning_array;
                bin_number = input_bin_number;
                binning = new TH1F(Form("%s binning", variable.Data()), "", bin_number, input_binning_array);

                AFBpos = new TH1F(Form("%spos", variable.Data()), Form("%spos", variable.Data()), bin_number, binning_array);
                AFBneg = new TH1F(Form("%sneg", variable.Data()), Form("%sneg", variable.Data()), bin_number, binning_array);
                AFBpos->Sumw2();
                AFBneg->Sumw2();
                mean = new TH1F(Form("%smean", variable.Data()), Form("%smean", variable.Data()), bin_number, 0, bin_number);
                meanN = new TH1F(Form("%smeanN", variable.Data()), Form("%smeanN", variable.Data()), bin_number, 0, bin_number);
        }

        void SetLimitsAngularBin(float input_phi_min, float input_phi_max, float input_theta_min, float input_theta_max)
        {
                phi_min = input_phi_min;
                phi_max = input_phi_max;
                theta_min = input_theta_min;
                theta_max = input_theta_max;
        }*/
};

#endif
