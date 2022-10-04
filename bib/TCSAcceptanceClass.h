#ifndef TCSAcceptanceclass
#define TCSAcceptanceclass

class Acceptance
{
public:
        int bin_t;
        int bin_Eg;
        int bin_M;
        int bin_Phi;
        int bin_Theta;
        TH2D **Acc;
        TFile *acceptanceFile;
        TH2D **AcceptancetM2forbinning;
        TH1D *AcceptanceEg;

        Acceptance(TString input_file, int input_bin_t, int input_bin_Eg, int input_bin_M, int input_bin_Phi, int input_bin_Theta)
        {
                acceptanceFile = new TFile(input_file);
                AcceptanceEg = (TH1D *)acceptanceFile->Get("AcceptanceEg");
                TH2D *AcceptancetM2forbinningHist = (TH2D *)acceptanceFile->Get("Egbin0");

                bin_t = input_bin_t;
                bin_Eg = input_bin_Eg;
                bin_M = input_bin_M;

                // Check binning t, Eg, M
                int nbBinsInT = AcceptancetM2forbinningHist->GetNbinsX();
                int nbBinsInM = AcceptancetM2forbinningHist->GetNbinsY();
                int nbBinsInEg = AcceptanceEg->GetNbinsX();

                if ((nbBinsInT != input_bin_t) || (nbBinsInM != input_bin_M) || (nbBinsInEg != input_bin_Eg))
                {
                        throw std::invalid_argument("Binning of the acceptance file does not match input binning (t,Eg,M)");
                }

                // Binning in t and M2
                AcceptancetM2forbinning = new TH2D *[nbBinsInEg];
                for (int Egbin = 0; Egbin < nbBinsInEg; Egbin++)
                {
                        AcceptancetM2forbinning[Egbin] = AcceptancetM2forbinningHist;
                        AcceptancetM2forbinning[Egbin]->SetName(Form("Egbin%d", Egbin));
                        AcceptancetM2forbinning[Egbin]->SetTitle(Form("Eg bin %d;-t;M2", Egbin));
                }

                // Theta/Phi plots for each t/Eg/M2 bin
                Acc = new TH2D *[nbBinsInT * nbBinsInM * nbBinsInEg];

                for (int Egbin = 0; Egbin < nbBinsInEg; Egbin++)
                {
                        for (int mass = 0; mass < nbBinsInM; mass++)
                        {
                                for (int tbinning = 0; tbinning < nbBinsInT; tbinning++)
                                {
                                        TString numt = TString::Itoa(tbinning, 10);
                                        TString numm = TString::Itoa(mass, 10);
                                        TString numE = TString::Itoa(Egbin, 10);
                                        TString name = numm + numt + numE;
                                        Acc[tbinning + nbBinsInT * mass + nbBinsInT * nbBinsInM * Egbin] = (TH2D *)acceptanceFile->Get(name);
                                }
                        }
                }

                // Check binning t, Eg, M
                int binInPhi = Acc[0]->GetNbinsX();
                int binInTheta = Acc[0]->GetNbinsY();

                bin_Phi = binInPhi;
                bin_Theta = binInTheta;

                if ((binInPhi != input_bin_Phi) || (binInTheta != input_bin_Theta))
                {
                        throw std::invalid_argument("Binning of the acceptance file does not match input binning (phi,theta)");
                }

                cout << "///////" << endl;
                cout << "Bin acceptance" << endl;
                cout << "///////" << endl;
                cout << " " << nbBinsInEg << " " << nbBinsInM << " " << nbBinsInT << " " << binInPhi << " " << binInTheta << endl;
                cout << "///////" << endl;
        }

        double GetValue(double t, double Epho, double qp2, double phi, double theta)
        {
                int binningEg = (AcceptanceEg->GetXaxis()->FindBin(Epho)) - 1;
                int binningt = (AcceptancetM2forbinning[binningEg]->GetXaxis()->FindBin(-t)) - 1;
                int binningMass = (AcceptancetM2forbinning[binningEg]->GetYaxis()->FindBin(qp2)) - 1;
                int binningphi = Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetXaxis()->FindBin(phi);
                int binningtheta = Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetYaxis()->FindBin(theta);

                double AccValue = Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetBinContent((Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetXaxis()->FindBin(phi)), (Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetYaxis()->FindBin(theta)));
                return AccValue;
        }

        double GetError(double t, double Epho, double qp2, double phi, double theta)
        {
                int binningEg = (AcceptanceEg->GetXaxis()->FindBin(Epho)) - 1;
                int binningt = (AcceptancetM2forbinning[binningEg]->GetXaxis()->FindBin(-t)) - 1;
                int binningMass = (AcceptancetM2forbinning[binningEg]->GetYaxis()->FindBin(qp2)) - 1;
                int binningphi = Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetXaxis()->FindBin(phi);
                int binningtheta = Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetYaxis()->FindBin(theta);

                double AccError = Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetBinError((Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetXaxis()->FindBin(phi)), (Acc[binningt + bin_t * binningMass + bin_t * bin_M * binningEg]->GetYaxis()->FindBin(theta)));
                return AccError;
        }

        int Get_t_bin(double t)
        { // to change if t binning become Eg-dependant
                return ((AcceptancetM2forbinning[0]->GetXaxis()->FindBin(-t)) - 1);
        }

        int Get_Eg_bin(double Epho)
        {
                return ((AcceptanceEg->GetXaxis()->FindBin(Epho)) - 1);
        }

        int Get_M_bin(double qp2)
        {
                // to change if M binning become Eg-dependant
                return ((AcceptancetM2forbinning[0]->GetYaxis()->FindBin(qp2)) - 1);
        }

        int Get_Phi_bin(double phi)
        {
                return (Acc[0]->GetXaxis()->FindBin(phi));
        }

        int Get_Theta_bin(double theta)
        {
                return (Acc[0]->GetYaxis()->FindBin(theta));
        }

        void Draw_Acc()
        {
               /* TCanvas *AcccanNewAcc0 = new TCanvas("AcccanNewAcc", "", 7000, 5000);
                AcccanNewAcc0->Divide(Acc_TCS.bin_t, Acc_TCS.bin_M);
                for (int mass = 0; mass < Acc_TCS.bin_M; mass++)
                {
                        for (int tbinning = 0; tbinning < Acc_TCS.bin_t; tbinning++)
                        {
                                int Egbin = 2;
                                AcccanNewAcc0->cd(1 + tbinning + Acc_TCS.bin_t * mass);
                                Acc_TCS.Acc[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin]->GetYaxis()->SetTitle("#theta (#circ)");
                                Acc_TCS.Acc[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin]->GetXaxis()->SetTitle("#phi (#circ)");
                                Acc_TCS.Acc[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin]->GetYaxis()->SetTitleOffset(0.9);
                                Acc_TCS.Acc[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin]->GetZaxis()->SetLabelSize(0.025);
                                Acc_TCS.Acc[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin]->Draw("colz");
                        }
                }
                AcccanNewAcc0->SaveAs("acccheck0.pdf");
                Acc_TCS.Acc[2 + Acc_TCS.bin_t * 1 + Acc_TCS.bin_t * Acc_TCS.bin_M * 2]->SaveAs("AccPlotNice.root");*/
        }

        void Draw_Error()
        {

              /*  TH1F *AcceptanceError = new TH1F("errorAcc", "; Error;nb bins", 100, 0, 1);
                TCanvas *AcceptanceErrorAcc = new TCanvas("AcceptanceErrorAcc", "", 6500, 5000);
                for (int Egbin = 0; Egbin < Acc_TCS.bin_Eg; Egbin++)
                {
                        for (int mass = 0; mass < Acc_TCS.bin_M; mass++)
                        {
                                for (int tbinning = 0; tbinning < Acc_TCS.bin_t; tbinning++)
                                {
                                        for (int phiBin = 1; phiBin < 37; phiBin++)
                                        {
                                                for (int thetabin = 1; thetabin < 14; thetabin++)
                                                {
                                                        double AccValue = Acc_TCS.Acc[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin]->GetBinContent(phiBin, thetabin);
                                                        double AccError = Acc_TCS.Acc[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin]->GetBinError(phiBin, thetabin);
                                                        if (AccValue > 0.05)
                                                                AcceptanceError->Fill(AccError / AccValue);
                                                }
                                        }
                                }
                        }
                }
                AcceptanceError->Draw();
                AcceptanceErrorAcc->SaveAs("AcceptanceErrorAcc.pdf");
                AcceptanceErrorAcc->SaveAs("AcceptanceErrorAcc.root");*/
        }
};

#endif
