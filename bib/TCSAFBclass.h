#ifndef TCSAFBclass
#define TCSAFBclass

class TCSAFB
{
public:
        TString variable;     // Name of the variable
        int bin_number;       // bin nummber
        TH1F *binning;        // binning Th1F
        float *binning_array; // Array with bin limit
        TH1F *AFBpos;         // Positive bin events
        TH1F *AFBneg;         // Negative bin events
        TH1F *mean;           // Store numerator of mean variable for each bin
        TH1F *meanN;          // Store denominator of mean variable for each bin

        float phi_min = -40.;
        float phi_max = 40.;

        float theta_min = 50.;
        float theta_max = 80.;

        float cut_acc = 0.05;
        float cut_acc_error = 0.5;

        TString outputfolder = "Plots"; // outputfolder for plots and graphs

        TCSAFB(TString input_variable, int input_bin_number, float *input_binning_array)
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
        }

        void SetCutAcc(float input_cut_acc, float input_cut_acc_error)
        {
                cut_acc = input_cut_acc;
                cut_acc_error = input_cut_acc_error;
        }

        void SetOutputFolder(TString input_outputfolder)
        {
                outputfolder = input_outputfolder;
        }

        void AddEvent(float value, float phi, float theta, float AccValue, float AccError, float weight, float CorrVolume1, float CorrVolume2)
        {

                int binning_var = (binning->GetXaxis()->FindBin(value)) - 1;
                if (AccValue > cut_acc && (AccError / AccValue) < cut_acc_error)
                {

                        if (phi > phi_min && phi < phi_max && theta < theta_max && theta > theta_min)
                        {
                                AFBpos->Fill(value, weight / (AccValue * CorrVolume1));
                                mean->Fill(binning_var, value * (weight / (AccValue * CorrVolume1)));
                                meanN->Fill(binning_var, (weight / (AccValue * CorrVolume1)));
                        }

                        if (theta < (180. - theta_min) && theta > (180. - theta_max) && (phi > (180. + phi_min) || phi < (-180. + phi_max)))
                        {
                                AFBneg->Fill(value, weight / (AccValue * CorrVolume2));
                                mean->Fill(binning_var, value * (weight / (AccValue * CorrVolume2)));
                                meanN->Fill(binning_var, (weight / (AccValue * CorrVolume2)));
                        }
                }
        }

        void Calculate()
        {

                TH1F *num = (TH1F *)AFBpos->Clone(Form("%spos", variable.Data()));
                TH1F *denom = (TH1F *)AFBpos->Clone(Form("%spos", variable.Data()));
                num->Sumw2();
                denom->Sumw2();

                num->Add(AFBneg, -1);
                denom->Add(AFBneg);
                num->Divide(denom);
                TGraphAsymmErrors *AFB_Graph = new TGraphAsymmErrors(bin_number);
                for (int i = 0; i < bin_number; i++)
                {
                        double bincenter = (mean->GetBinContent(i + 1)) / (meanN->GetBinContent(i + 1));
                        AFB_Graph->SetPoint(i, bincenter, (num->GetBinContent(i + 1)));
                        AFB_Graph->SetPointError(i, (bincenter - binning_array[i]), (binning_array[i + 1] - bincenter), (num->GetBinError(i + 1)), (num->GetBinError(i + 1)));
                }
                AFB_Graph->SetMarkerSize(5);
                AFB_Graph->SetMarkerColor(kRed);
                AFB_Graph->SetLineColor(kRed);
                AFB_Graph->GetXaxis()->SetRangeUser(1.5, 3.);
                AFB_Graph->SetTitle(Form(";%s;A_{FB}", variable.Data()));

                TCanvas *canAFB = new TCanvas(Form("canAFB%s", variable.Data()), Form("canAFB%s", variable.Data()), 5000, 1000);
                canAFB->Divide(4, 1);
                canAFB->cd(1);
                AFBpos->Draw();
                canAFB->cd(2);
                AFBneg->Draw();
                canAFB->cd(3);
                num->Draw();
                canAFB->SaveAs(Form("%s/canAFB%s.pdf", outputfolder.Data(), variable.Data()));
        }
};

#endif
