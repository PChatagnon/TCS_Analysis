#ifndef TCSPlotClass
#define TCSPlotClass

const std::vector<vector<TString>> ListPlot1D{

    // Format : { "name", "title;tile_x;title_y" , bin_x , min_x , max_x }
    {"evt_count", "evt_count", "4", "0", "4"},
    {"efficiency", "efficiency", "3", "0", "3"},
    {"MMass", "Missing mass beam", "100", "-20", "20"},
    {"Q2", "Virtuality incoming photon (Q2)", "100", "-20", "20"},
    {"EM1", ";M;events", "200", "0.0", "4.0"},

    {"Cons. Energy", "Cons. Energy", "100", "-1", "1"},
    {"T rec", "t REC", "50", "-0.2", "2"},

    {"RECpid", "Histogram of rec pid", "100", "-15", "2500"},
    {"MMass", "Missing mass beam", "100", "-20", "20"},
    {"Q2", "Virtuality incoming photon (Q2)", "100", "-20", "20"},
    {"M2", "Invariant Mass of the lepton pair squared", "120", "-1", "5"},
    {"M", "Invariant Mass of the lepton pair", "200", "0", "4"},
    {"xihist1", "Eg(4-10 GeV),t[0.15-0.8];#xi;Q2", "100", "0.05", "0.35"},

    {"ET rec", "t REC", "50", "-0.2", "2"},

    {"EMMass", "Missing mass beam", "100", "-1", "1"},
    {"EQ2", "Virtuality incoming photon (Q2)", "80", "-0.1", "0.3"},
    {"EQ21", "Virtuality incoming photon (Q2)", "80", "-0.1", "0.3"},
    {"EM", "Invariant Mass of the lepton pair", "120", "-1", "10"},
    {"EM1", "", "200", "0.0", "4.0"},
    {"EM1AfterFid", "Invariant Mass of the lepton pair AfterFid", "200", "0", "4"},
    {"EM1prim", "Invariant Mass of the lepton pair", "200", "0", "3"},
    {"EM2", "Invariant Mass of the lepton pair with scattered electron", "120", "0", "3"},

    {"ECons. Energy", "Cons. Energy", "100", "-1", "1"},
    {"ECons. EnergyEvents", "Cons. Energy", "100", "-1", "1"},

    {"FinalEventst", "", "50", "0", "1"},
    {"FinalEventsEg", "", "50", "0", "11"},
    {"FinalEventsM", "", "160", "0.", "4"},
    {"FinalEventsPhi", "", "50", "-180", "180"},
    {"FinalEventsP", "", "50", "0", "1.2"},
    {"FinalEventsTheta", "", "50", "5", "60"},
    {"EvertexPair", "Vertex diff pair", "50", "-10", "10"},
    {"EvertexElecP", "Vertex Elec p", "50", "-10", "10"},
    {"EvertexPosP", "Vertex Posi p", "50", "-10", "10"},
    {"vertexTimePair", "vertex time pair", "150", "-2", "2"},
    {"vertexTimeP", "vertex time proton", "150", "-2", "2"},
    {"Pproton", "Pproton", "100", "0.", "100"},
    {"EvertexTimePair", "Vertex time pair", "150", "-2", "2"},
    {"EvertexTimeP", "vertex time proton", "100", "-0.5", "0.5"},

    {"vertexPair", "Vertex diff pair", "50", "-10", "10"},
    {"vertexElecP", "Vertex Elec p", "50", "-10", "10"},
    {"vertexPosP", "Vertex Posi p", "50", "-10", "10"},
    {"vertexElec", "Vertex elec", "50", "-10", "10"},
    {"vertexPosi", "Vertex posi ", "50", "-10", "10"},
    {"vertexProt", "Vertex proton", "50", "-10", "10"},
    {"Chi2ElectronAvant", "Electron (Before exclu.} #chi^{2}", "100", "-6", "6"},
    {"Chi2PositronAvant", "Positron (Before exclu.} #chi^{2}", "100", "-6", "6"},
    {"Chi2ProtonAvant", "Proton CD (Before exclu.} #chi^{2}", "100", " -10", "10"},
    {"Chi2ProtonFDAvant", "Proton FD (Before exclu.} #chi^{2}", "100", "-10", "10"},

    {"Chi2Electron", "Electron; #chi^{2}", "100", "-6", "6"},
    {"Chi2Positron", "Positron; #chi^{2}", "100", "-6", "6"},
    {"Chi2Proton", "Proton CD; #chi^{2}", "100", "-10", "10"},
    {"Chi2ProtonFD", "Proton FD; #chi^{2}", "100", "-10", "10"},
    {"MLPplot", "", "100", "-0.05", "1.05"},
    {"PositronPt", "PositronPt", "100", "0.", "4."},
    {"tCheck", "", "100", "-1", "1"},
    {"VertexOtherPart", "VertexOtherPart", "500", "-30", "30"},
    {"Chi2OtherPart", "Chi2OtherPart", "500", "-40", "40"},

    {"tCheckEvents", "", "100", "-1", "1"},

    {"tCheckBefore", "", "100", "-1", "1"},
    {"Cherenkov Electron", "Cherenkov Electron", "40", "0", "40"},
    {"Cherenkov Positron", "Cherenkov Positron", "40", "0", "40"},
    {"ECherenkov Electron", "Cherenkov Electron", "40", "0", "40"},
    {"ECherenkov Positron", "Cherenkov Positron", "40", "0", "40"}};

const std::vector<vector<TString>> ListPlot2D{

    // Format : {"name", "title;tile_x;title_y",bin_x,min_x,max_x,bin_y,min_y,max_y}
    {"MvsE", "M vs Egamma", "120", "2", "11", "120", "0", "3"},
    {"ChePhi", "ChePhi", "200", "-180", "180", "200", "0", "40"},
    {"MvsPElec", "M vs PElec", "120", "0", "11", "120", "0", "3"},
    {"thist", "Mass[1.5-2],Eg(4-10 GeV),t[0.15-0.8];#xi;-t", "100", "0.05", "0.35", "100", "0.15", "0.8"},
    {"thist1", "Mass[2-3],Eg(4-10 GeV),t[0.15-0.8];#xi;-t", "100", "0.05", "0.35", "100", "0.15", "0.8"},
    {"Before Cuts", " Before Cuts", "200", "-400", "400", "200", "-400", "400"},
    {"Before Cuts Positron", " Before Cuts", "200", "-400", "400", "200", "-400", "400"},

    {"MvsE", "M vs Egamma", "120", "2", "11", "120", "0", "3"},
    {"MvsPElec", "M vs PElec", "120", "0", "11", "120", "0", "3"},
    {"MvsPPosi", "M vs PPosi", "120", "0", "11", "120", "0", "3"},
    {"MMvsE", "MM vs Egamma", "100", "2", "11", "100", "-5", "10"},
    {"PMiss", "Direction of the missing particle", "50", "-0.5", "0.5", "50", "-0.5", "0.5"},
    {"MMPt", "T. Momentum Fraction Vs MM beam", "100", "-1", "1", "50", "0", "0.1"},
    {"xihist", "Eg(4-10 GeV),t[0.15-0.8];#xi;Q2", "100", "0.05", "0.35", "100", "2.", "9."},

    {"MMQ2", " Q2 VS Missing mass beam", "50", "-3", "3", "50", "-3", "3"},
    {"RECq2vst", "Q2 vs -t (REC Particles)", "50", "4", "9", "30", "0", "1"},
    {"PhiVSPelectron", "P vs Phi (Electron)", "100", "-180", "180", "100", "0", "10"},
    {"PhiVSPPositron", "P vs Phi (Positron)", "100", "-180", "180", "100", "0", "10"},
    {"PhiVSPProton", "P vs Phi (Proton)", "100", "-180", "180", "100", "0", "5"},
    {"ThetaVSPelectron", "Theta vs P (Electron)", "100", "0", "10", "100", "0", "50"},
    {"ThetaVSPPositron", "Theta vs P (Positron)", "100", "0", "10", "100", "0", "50"},
    {"MvsPhiProton", "M vs PhiProton", "120", "-190", "190", "120", "0", "3"},
    {"MvsPhiElectron", "M vs PhiElec", "120", "-190", "190", "120", "0", "3"},
    {"MvsPhiPositron", "M vs PhiPosi", "120", "-190", "190", "120", "0", "3"},
    {"MvsThetaElectrogdgdgn", "M vs ThetaElec", "100", "0", "50", "120", "0", "3"},
    {"MvsThetaElectron", "M vs ThetaElec", "100", "0", "0", "120", "0", "3"},
    {"MvsThetaPositron", "M vs ThetaPosi", "100", "0", "50", "120", "0", "3"},
    {"MvsThetaProton", "M vs ThetaProt", "100", "10", "80", "120", "0", "3"},
    {"MvsVElectron", "M vs Vx diff pair", "100", "-10", "10", "120", "0", "3"},
    {"MvsVPositron", "M vs Vy diff pair", "100", "-10", "10", "120", "0", "3"},
    {"MvsVProton", "M vs Vz diff pair", "100", "-10", "10", "120", "0", "3"},
    {"MvsCheElectron", "M vs CheElec", "100", "0", "40", "120", "0", "3"},
    {"MvsChePositron", "M vs ChePosi", "100", "0", "40", "120", "0", "3"},
    {"MvsCheProton", "M vs CheProt", "100", "0", "40", "120", "0", "3"},
    {"MvsCheElectron", "ChePosi vs CheElec", "100", "0", "40", "100", "0", "40"},
    {"MvsCheElectrondgdgdd", "ChePosi vs CheElec", "100", "0", "40", "100", "0", "40"},
    {"MvsCheElectronaaaa", "ChePosi vs CheElec", "100", "0", "40", "100", "0", "40"},
    {"PhiVSThetaelectron", "Theta vs Phi (Electron)", "100", "-180", "180", "100", "-0.1", "50"},
    {"PhiVSThetaPositron", "Theta vs Phi (Positron)", "100", "-180", "180", "100", "-0.1", "50"},
    {"ThetaVSThetaProton", "Theta vs Phi (Proton)", "100", "-180", "180", "100", "-0.1", "70"},
    {"ThetaVSPProton", "Theta vs P (Proton)", "100", "0", "5", "100", "0", "70"},

    {"EM1SD", "Sector Diff vs Invariant Mass of the lepton pair", "120", "0", "3", "6", "-1", "5"},
    {"EMvsE", "M vs Egamma", "100", "2", "11", "100", "0", "2.5"},
    {"EMvsPElec", "M vs P Electron", "100", "0", "11", "100", "0", "2.5"},
    {"EMvsPPosi", "M vs P Positron", "100", "0", "11", "100", "0", "2.5"},
    {"EMMvsE", "MM vs Egamma", "100", "2", "11", "100", "-5", "10"},
    {"EMvsPProton", "M vs P Proton", "100", "0", "2", "100", "0", "2.5"},
    {"EMvsRun", "M vs RUN", "65", "4013", "4078", "30", "0", " 2.5"},
    {"EMvsPhiProton", "M vs PhiProton", "120", "-190", "190", "120", "0", "3"},
    {"EPMiss", "Direction of the missing particle", "50", "-0.5", "0.5", "50", "-0.5", "0.5"},
    {"EMMPt", "Pt Vs MM", "100", "-5", "5", "50", "0", "1"},

    {"EMMQ2", "MM Q2", "50", "-15", "15", "50", "-6", "6"},
    {"EPhiVSPelectron", "P vs Phi (Electron)", "100", "-180", "180", "100", "0", "10"},
    {"EPhiVSPPositron", "P vs Phi (Positron)", "100", "-180", "180", "100", "0", "10"},
    {"EPhiVSPProton", "P vs Phi (Proton)", "100", "-180", "180", "100", "0", "5"},
    {"EThetaVSPelectron", "Theta vs P (Electron)", "100", "0", "10", "100", "0", "50"},
    {"EThetaVSPPositron", "Theta vs P (Positron)", "100", "0", "10", "100", "0", "50"},
    {"EPhiVSThetaelectron", "Theta vs Phi (Electron)", "100", "-180", "180", "100", "-0.1", "50"},
    {"EPhiVSThetaPositron", "Theta vs Phi (Positron)", "100", "-180", "180", "100", "-0.1", "50"},
    {"EThetaVSThetaProton", "Theta vs Phi (Proton)", "100", "-180", "180", "100", "-0.1", "70"},
    {"EThetaVSPProton", "Theta vs P (Proton)", "100", "0", "5", "50", "0", "70"},

    {"vertexProtTheta", "Vertex proton vs Theta", "100", "0", "70", "50", "-10", "10"},
    {"vertexProtThetaCD", "Vertex proton vs Theta CD", "100", "0", "70", "50", "-10", "10"},
    {"vertexProtThetaFD", "Vertex proton vs Theta FD", "100", "0", "70", "50", "-10", "10"},

    {"vertexTimePP", "vertex time proton vs p proton", "100", "0", "3", "100", "-2", "2"},

    {"vertexTimePairMass", "", "100", "0", "3", "150", "-2", "2"},

    {"EvertexTimePP", "vertex time proton vs p proton CD", "60", "0.25", "1", "100", "-0.5", "0.5"},
    {"EvertexTimePPFD", "vertex time proton vs p proton FD", "60", "0.25", "1", "100", "-0.5", "0.5"},

    {"EThetaVSPelectron", "Theta vs P (Electron)", "100", "0", "10", "100", "0", "50"},
    {"EThetaVSPPositron", "Theta vs P (Positron)", "100", "0", "10", "100", "0", "50"},
    {"EThetaVSPProton", "Theta vs P (Proton)", "50", "0", "1.5", "50", "0", "70"},
    {"PPosiPElec", "All events;P (e+},P (e-)", "50", "0", "10", "50", "0", "10"},
    {"PPosiPElec1", "M>1.5 GeV;P (e+},P (e-)", "50", "0", "10", "50", "0", "10"},

    {"SFelectron", "SF vs P (Electron)", "80", "0", "8", "100", "0.", "0.35"},
    {"SFpositron", "SF vs P (Positron)", "80", "0", "8", "100", "0.", "0.35"},
    {"corrSFelectron", "corrSF vs P (Electron)", "80", "0", "8", "100", "0.1", "0.35"},
    {"corrSFpositron", "corrSF vs P (Positron)", "80", "0", "8", "100", "0.1", "0.35"},
    {"ECelectron", "ECout vs ECin (Electron)", "50", "0", "0.7", "50", "0", "0.7"},
    {"ECpositron", "ECout vs ECin (Positron)", "50", "0", "0.7", "50", "0", "0.7"},
    {"EECelectron", "ECout vs ECin (Electron)", "50", "0", "0.7", "50", "0", "0.7"},
    {"EECpositron", "ECout vs ECin (Positron)", "50", "0", "0.7", "50", "0", "0.7"},

    {"Beta P Proton", " Beta P Proton ", "100", "0", "2.5", "50", "0", "1.2"},
    {"EBeta P Proton", " Beta P Proton ", "100", "0", "2.5", "50", "0", "1.2"},

    {" SF Electron vs U", " SF Electron vs U", "150", "0", "450", "50", "0.", "0.35"},
    {" SF Electron vs V", " SF Electron vs V", "150", "0", "450", "50", "0.", "0.35"},
    {" SF Electron vs W", " SF Electron vs W", "150", "0", "450", "50", "0.", "0.35"},
    {" SF Positron vs U", " SF Positron vs U", "150", "0", "450", "50", "0.", "0.35"},
    {" SF Positron vs V", " SF Positron vs V", "150", "0", "450", "50", "0.", "0.35"},
    {" SF Positron vs W", " SF Positron vs W", "150", "0", "450", "50", "0.", "0.35"},

    {" FiducialPCAL", "FiducialPCAL", "500", "-500", "500", "500", "-500", "500"},
    {" EFiducialPCAL", "FiducialPCAL", "500", "-500", "500", "500", "-500", "500"},

    {"EgVSPeletron", "EgVSPeletron", "100", "0", "10", "100", "3", "10"},
    {"tVSPeletron", "tVSPeletron", "100", "0", "10", "100", "0.0", "1."},

    {"SFMass", "SFMassElec", "100", "0", "3", "100", "0.", "0.3"},
    {"SFMass", "SFMassPosi", "100", "0", "3", "100", "0.", "0.3"},

    {"EgVST", "EgVST", "80", "3", "11", "80", "0.05", "1.1"},
    {"MVST", "MVST", "80", "1.4", "3.1", "80", "0.05", "1.1"},
    {"MVST", "MVST", "80", "1.4", "3.1", "80", "2", "11"},

    {"PhiMissingParticle", ";phi;Pt/P", "100", "-180.", "180", "100", "0", "0.7"},
    {"SFPCALECALPosi", "SFPCALECALPosi", "100", "0", "0.3", "100", "0", "0.3"},
    {"SFPCALECALPosia", "SFPCALECALPosia", "100", "0", "0.3", "100", "0", "0.3"},

    {"SFPCALECALElec", "SFPCALECALElec", "100", "0", "0.3", "100", "0.", "0.3"},
    {"SFPCALECALEleca", "SFPCALECALEleca", "100", "0", "0.3", "100", "0.", "0.3"},

    {"SFPpcalzeroPosi", "SFPpcalzeroPosi", "100", "0", "11", "100", "0.", "0.3"},
    {"SFPpcalzeroElec", "SFPpcalzeroElec", "100", "0", "11", "100", "0", "0.3"},

    {"tmass", ";t;mass", "100", "0.15", "0.8", "100", "2", "9."},
    {"PolarizationTransfer", "PolarizationTransfer", "100", "0", "11", "100", "0.", "1."}};

class TCSPlots
{
public:
        TString outputFolder = "Plots";
        std::unordered_map<string, TH1F *> Plot1D;
        std::unordered_map<string, TH2F *> Plot2D;

        void SetOutputFolder(TString input_outputFolder)
        {
                outputFolder = input_outputFolder;
        }

        void SavePlots(TString format)
        {
                for (auto plot : Plot1D)
                {
                        TCanvas *can = new TCanvas("", "can", 1500, 1000);
                        can->cd();
                        plot.second->Draw();
                        can->SaveAs(outputFolder + (plot.second->GetName()) + "." + format);
                }

                for (auto plot : Plot2D)
                {
                        TCanvas *can = new TCanvas("", "can", 1500, 1000);
                        can->cd();
                        plot.second->Draw("col");
                        can->SaveAs(outputFolder + (plot.second->GetName()) + "." + format);
                }
        }

        void Initialize_1D()
        {
                for (int i = 0; i < ListPlot1D.size(); i++)
                {
                        TString name = ListPlot1D[i][0];
                        TString title = ListPlot1D[i][1];
                        TString bin_x = ListPlot1D[i][2];
                        TString min_x = ListPlot1D[i][3];
                        TString max_x = ListPlot1D[i][4];

                        Plot1D[(string)name.Data()] = new TH1F(name, title, stoi((string)bin_x.Data()), stof((string)min_x.Data()), stof((string)max_x.Data()));
                }
        }

        void Initialize_2D()
        {
                for (int i = 0; i < ListPlot2D.size(); i++)
                {
                        TString name = ListPlot2D[i][0];
                        TString title = ListPlot2D[i][1];
                        TString bin_x = ListPlot2D[i][2];
                        TString min_x = ListPlot2D[i][3];
                        TString max_x = ListPlot2D[i][4];
                        TString bin_y = ListPlot2D[i][5];
                        TString min_y = ListPlot2D[i][6];
                        TString max_y = ListPlot2D[i][7];

                        Plot2D[(string)name.Data()] = new TH2F(name, title, stoi((string)bin_x.Data()), stof((string)min_x.Data()), stof((string)max_x.Data()), stoi((string)bin_y.Data()), stof((string)min_y.Data()), stof((string)max_y.Data()));
                }
        }

        void Fill_1D(TString name, double value, double weight)
        {
                Plot1D[(string)name.Data()]->Fill(value, weight);
        }

        void Fill_2D(TString name, double value_x, double value_y, double weight)
        {
                Plot2D[(string)name.Data()]->Fill(value_x, value_y, weight);
        }

        void Add_Hist_1D(TString name, TString title, int bin, float min_x, float max_x)
        {
                Plot1D[(string)name.Data()] = new TH1F(name, title, bin, min_x, max_x);
        }

        void Add_Hist_2D(TString name, TString title, int bin_x, float min_x, float max_x, int bin_y, float min_y, float max_y)
        {
                Plot2D[(string)name.Data()] = new TH2F(name, title, bin_x, min_x, max_x, bin_y, min_y, max_y);
        }

        /*   void Draw_Hist_1D(TString name)
           {
                   TCanvas *can = new TCanvas("canChi2", "canChi2", 4000, 2000);
           Chi2ElectronAvant->Draw();
           can->SaveAs("canChi2.pdf");
           can->SaveAs("canChi2.root");
           Chi2ProtonAvant->SaveAs("Chi2ProtonAvant.root");
           }

           void Draw_Hist_2D(TString name)
           {
                    TCanvas *can = new TCanvas("canChi2", "canChi2", 4000, 2000);
           Chi2ElectronAvant->Draw();
           can->SaveAs("canChi2.pdf");
           can->SaveAs("canChi2.root");
           Chi2ProtonAvant->SaveAs("Chi2ProtonAvant.root");
           }
   */
        void Draw_All_1D()
        {
        }

        void Draw_All_2D()
        {
        }
};

#endif
