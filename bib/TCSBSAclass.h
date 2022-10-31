#ifndef TCSBSAclass
#define TCSBSAclass

class TCSBSA
{
public:
        TString variable;     // Name of the variable
        int bin_number;       // bin nummber
        TH1F *binning;        // binning Th1F
        float *binning_array; // Array with bin limit
        TH1F **BSApos;        // Positive helicity events
        TH1F **BSAneg;        // Negative helicity events
        TH1F *meanBSA;        // Store numerator of mean variable for each bin
        TH1F *meanBSAN;       // Store denominator of mean variable for each bin

        TString outputfolder = "Plots"; // outputfolder for plots and graphs
        int bin_fit;                    // number of phi bin for the fit

        float cut_acc = 0.05;
        float cut_acc_error = 0.5;

        TCSBSA(TString input_variable, int input_bin_number, float *input_binning_array, int input_bin_fit)
        {
                variable = input_variable;
                bin_number = input_bin_number;
                binning = new TH1F(Form("%s binning", variable.Data()), "", bin_number, input_binning_array);
                bin_fit = input_bin_fit;
                binning_array = input_binning_array;

                BSApos = new TH1F *[bin_number];
                BSAneg = new TH1F *[bin_number];
                for (int bin = 0; bin < bin_number; bin++)
                {
                        BSApos[bin] = new TH1F(Form("%s bin %d pos", variable.Data(), bin), Form("%s bin %d pos", variable.Data(), bin), bin_fit, 0, 360);
                        BSAneg[bin] = new TH1F(Form("%s bin %d neg", variable.Data(), bin), Form("%s bin %d neg", variable.Data(), bin), bin_fit, 0, 360);
                        BSApos[bin]->Sumw2();
                        BSAneg[bin]->Sumw2();
                }

                meanBSA = new TH1F("meanBSA", "meanBSA", bin_number, 0, bin_number);
                meanBSAN = new TH1F("meanBSAN", "meanBSAN", bin_number, 0, bin_number);
        }

        void setOutputFolder(TString input_outputfolder)
        {
                outputfolder = input_outputfolder;
        }

         void SetCutAcc(float input_cut_acc, float input_cut_acc_error)
        {
                cut_acc = input_cut_acc;
                cut_acc_error = input_cut_acc_error;
        }

        void AddEvent(float value, float phi, int polarization, float polaT, float AccValue, float AccError, float weight)
        {
                int binning_var = (binning->GetXaxis()->FindBin(value)) - 1;
                // cout<<"Binning var "<<binning_var<<" value "<<value<<" "<<variable<<endl;
                if (phi < 0.0)
                        phi = phi + 360.;


                //cout<<value<<" "<< phi<<" "<< polarization<<" "<< polaT<<" "<< AccValue<<" "<< AccError<<" "<< weight<<endl;
                if (AccValue > cut_acc && (AccError / AccValue) < cut_acc_error)
                {
                        if (polarization == 1)
                                BSApos[binning_var]->Fill(phi, weight / (AccValue * polaT));
                        if (polarization == -1)
                                BSAneg[binning_var]->Fill(phi, weight / (AccValue * polaT));

                        meanBSA->Fill(binning_var, value * (weight / AccValue));
                        meanBSAN->Fill(binning_var, (weight / AccValue));
                }
        }

        void Calculate(float factorPola)
        {
                TF1 *f2 = new TF1("f2", "[0]*sin(x*[1]/[2])", 0, 360);
                f2->FixParameter(1, 3.14159264);
                f2->FixParameter(2, 180);
                f2->SetParNames("A_{#odot U}", "p1", "p2");

                TGraphAsymmErrors *BSAversus = new TGraphAsymmErrors(bin_number);
                TCanvas *canpola = new TCanvas("canpola", "canpola", 4000, 4750);
                canpola->Divide(2, 3);
                for (int bin = 0; bin < bin_number; bin++)
                {

                        TH1F *numpolap = (TH1F *)BSApos[bin]->Clone(Form("%s bin %d pos", variable.Data(), bin));
                        TH1F *denom1polap = (TH1F *)BSApos[bin]->Clone(Form("%s bin %d pos", variable.Data(), bin));
                        numpolap->Sumw2();
                        denom1polap->Sumw2();
                        numpolap->SetTitle(Form("Bin %d;#phi;BSA", bin));
                        numpolap->Add(BSAneg[bin], -1);
                        denom1polap->Add(BSAneg[bin]);
                        numpolap->Divide(denom1polap);
                        numpolap->Scale(1. / factorPola);
                        numpolap->GetXaxis()->SetTitle("#phi");
                        gStyle->SetOptFit(1);
                        gStyle->SetOptStat(0);
                        canpola->cd(bin + 1);
                        numpolap->Fit("f2");
                        numpolap->SaveAs(Form("%s/Fit%s%d.root", outputfolder.Data(), variable.Data(), bin));

                        double fitresult = numpolap->GetFunction("f2")->GetParameter(0);
                        double fiterror = numpolap->GetFunction("f2")->GetParError(0);

                        double bincenter = (meanBSA->GetBinContent(bin + 1)) / (meanBSAN->GetBinContent(bin + 1));
                        BSAversus->SetPoint(bin, bincenter, fitresult);
                        BSAversus->SetPointError(bin, (bincenter - binning_array[bin]), (binning_array[bin + 1] - bincenter), fiterror, fiterror);
                }
                canpola->cd(5);
                BSAversus->SetFillColor(kRed);
                BSAversus->SetFillStyle(3005);
                BSAversus->SetMarkerSize(3);
                BSAversus->SetMarkerColor(kRed);
                BSAversus->SetLineColor(kRed);

                BSAversus->SetTitle(Form(";%s;BSA", variable.Data()));
                BSAversus->Draw("A P");
                canpola->SaveAs(outputfolder + "/canpola" + variable + ".pdf");
                canpola->SaveAs(outputfolder + "/canpola" + variable + ".root");
                BSAversus->SaveAs(outputfolder + "/BSAversus" + variable + ".root");
        }
};

#endif
