#ifndef TCSRRatioclass
#define TCSRRatioclass

class TCSRRatio
{
public:
        TString variable;     // Name of the variable
        int bin_number;       // bin nummber
        TH1F *binning;        // binning Th1F
        float *binning_array; // Array with bin limit
        TH1F **RratioNum;
        TH1F **RratioDenom;
        TH1F *mean;  // Store numerator of mean variable for each bin
        TH1F *meanN; // Store denominator of mean variable for each bin

        float cut_acc = 0.05;
        float cut_acc_error = 0.5;

        TString outputfolder = "Plots"; // outputfolder for plots and graphs

        TCSRRatio(TString input_variable, int input_bin_number, float *input_binning_array)
        {
                variable = input_variable;
                binning_array = input_binning_array;
                bin_number = input_bin_number;
                binning = new TH1F(Form("%sbinning_R", variable.Data()), "", bin_number, input_binning_array);

                RratioNum = new TH1F*[bin_number];
                RratioDenom = new TH1F*[bin_number];

                for (int bin = 0; bin < bin_number; bin++)
                {
                        RratioNum[bin] = new TH1F(Form("%s bin %d num", variable.Data(), bin), Form("%s bin %d num", variable.Data(), bin), 10, -180, 180);
                        RratioDenom[bin] = new TH1F(Form("%s bin %d denom", variable.Data(), bin), Form("%s bin %d denom", variable.Data(), bin), 10, -180, 180);
                        RratioNum[bin]->Sumw2();
                        RratioDenom[bin]->Sumw2();
                }

                mean = new TH1F(Form("%smean_R", variable.Data()), Form("%smean_R", variable.Data()), bin_number, 0, bin_number);
                meanN = new TH1F(Form("%smeanN_R", variable.Data()), Form("%smeanN_R", variable.Data()), bin_number, 0, bin_number);
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

        void AddEvent(float value, float phi, float theta, float L, float L0, float AccValue, float AccError, float weight)
        {
                int binning_var = (binning->GetXaxis()->FindBin(value)) - 1;
                if (AccValue > cut_acc && (AccError / AccValue) < cut_acc_error)
                {
                        if (theta > 45. && theta < 135.)
                        {

                                RratioNum[binning_var]->Fill(phi, (L / L0) * (weight / AccValue));
                                RratioDenom[binning_var]->Fill(phi, (L / L0) * (weight / AccValue));
                                mean->Fill(binning_var, value * (weight / AccValue));
                                meanN->Fill(binning_var, weight / AccValue);
                        }
                }
        }

        void Calculate()
        {

                TF1 *f3 = new TF1("f2", "cos(x*[1]/[2])", -180, 180);
                f3->FixParameter(1, 3.14159264);
                f3->FixParameter(2, 180);
                TGraphAsymmErrors *Rratio = new TGraphAsymmErrors(bin_number);
                TCanvas *canRatio = new TCanvas("canRatio", "canRatio", (bin_number*2000), 4000);
                canRatio->Divide(bin_number, 3);
                for (int bin = 0; bin < bin_number; bin++)
                {

                        TH1F *numR = (TH1F *)RratioNum[bin]->Clone(Form("%s bin %d num_R", variable.Data(), bin));
                        TH1F *denomR = (TH1F *)RratioDenom[bin]->Clone(Form("%s bin %d denom_R", variable.Data(), bin));

                        numR->Multiply(f3, 1.);

                        canRatio->cd(bin + 1);
                        numR->Draw();
                        canRatio->cd(bin + 1 + bin_number);
                        RratioNum[bin]->Draw();

                        double RC = numR->Integral();
                        double R0 = denomR->Integral();
                        double R = RC / R0;

                        double errorY = smear(*RratioNum[bin], bin);
                        double bincenter = (mean->GetBinContent(bin + 1)) / (meanN->GetBinContent(bin + 1));
                        Rratio->SetPoint(bin, bincenter, R);
                        Rratio->SetPointError(bin, (bincenter - binning_array[bin]), (binning_array[bin + 1] - bincenter), errorY, errorY);
                }
                canRatio->cd((bin_number*2)+1);
                Rratio->SetFillColor(kRed);
                Rratio->SetFillStyle(3005);
                Rratio->SetMarkerSize(3);
                Rratio->SetMarkerColor(kRed);
                Rratio->SetLineColor(kRed);
                Rratio->SetTitle(Form(";%s;R'", variable.Data()));
                Rratio->Draw("A P");
                canRatio->SaveAs(outputfolder + "/canRatio" + variable + ".pdf");
                Rratio->SaveAs(outputfolder + "/Ratio" + variable + ".root");
        }
};

#endif
