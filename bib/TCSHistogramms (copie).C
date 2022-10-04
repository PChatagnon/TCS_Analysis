#include "TCSHistogramms.h"

FinalEventst = new TH1F("FinalEventst", "", 50, 0, 1);
FinalEventsEg = new TH1F("FinalEventsEg", "", 50, 0, 11);
FinalEventsM = new TH1F("FinalEventsM", "", 160, 0., 4);
FinalEventsPhi = new TH1F("FinalEventsPhi", "", 50, -180, 180);
FinalEventsP = new TH1F("FinalEventsP", "", 50, 0, 1.2);
FinalEventsTheta = new TH1F("FinalEventsTheta", "", 50, 5, 60);
xihist = new TH2F("xihist", "Eg(4-10 GeV),t[0.15-0.8];#xi;Q2", 100, 0.05, 0.35, 100, 2., 9.);
xihist1 = new TH1F("xihist1", "Eg(4-10 GeV),t[0.15-0.8];#xi;Q2", 100, 0.05, 0.35);
thist = new TH2F("thist", "Mass[1.5-2],Eg(4-10 GeV),t[0.15-0.8];#xi;-t", 100, 0.05, 0.35, 100, 0.15, 0.8);
thist1 = new TH2F("thist1", "Mass[2-3],Eg(4-10 GeV),t[0.15-0.8];#xi;-t", 100, 0.05, 0.35, 100, 0.15, 0.8);
BeforeCuts = new TH2F("Before Cuts", " Before Cuts", 200, -400, 400, 200, -400, 400);
BeforeCutsPositron = new TH2F("Before Cuts Positron", " Before Cuts", 200, -400, 400, 200, -400, 400);
histRECPid = new TH1F("RECpid", "Histogram of rec pid", 100, -15, 2500);
histMMass = new TH1F("MMass", "Missing mass beam", 100, -20, 20);
histQ2 = new TH1F("Q2", "Virtuality incoming photon (Q2)", 100, -20, 20);
histQP2 = new TH1F("M2", "Invariant Mass of the lepton pair squared", 120, -1, 5);
histM = new TH1F("M", "Invariant Mass of the lepton pair", 200, 0, 4);
histME = new TH2D("MvsE", "M vs Egamma", 120, 2, 11, 120, 0, 3);
histMPElec = new TH2D("MvsPElec", "M vs PElec", 120, 0, 11, 120, 0, 3);
histMPPosi = new TH2D("MvsPPosi", "M vs PPosi", 120, 0, 11, 120, 0, 3);
histMME = new TH2D("MMvsE", "MM vs Egamma", 100, 2, 11, 100, -5, 10);
histPMiss = new TH2D("PMiss", "Direction of the missing particle", 50, -0.5, 0.5, 50, -0.5, 0.5);
histT = new TH1F("T rec", "t REC", 50, -0.2, 2);
histMMPt = new TH2D("MMPt", "T. Momentum Fraction Vs MM beam", 100, -1, 1, 50, 0, 0.1); //-1,1,50,0,0.2);//
EnergyConservation = new TH1F("Cons. Energy", "Cons. Energy", 100, -1, 1);
MMQ2 = new TH2D("MMQ2", " Q2 VS Missing mass beam", 50, -3, 3, 50, -3, 3);
histQ2t = new TH2D("RECq2vst", "Q2 vs -t (REC Particles)", 50, 4, 9, 30, 0, 1);
PhiVSPElectron = new TH2D("PhiVSPelectron", "P vs Phi (Electron)", 100, -180, 180, 100, 0, 10);
PhiVSPPositron = new TH2D("PhiVSPPositron", "P vs Phi (Positron)", 100, -180, 180, 100, 0, 10);
PhiVSPProton = new TH2D("PhiVSPProton", "P vs Phi (Proton)", 100, -180, 180, 100, 0, 5);
ThetaVSPElectron = new TH2D("ThetaVSPelectron", "Theta vs P (Electron)", 100, 0, 10, 100, 0, 50);
ThetaVSPPositron = new TH2D("ThetaVSPPositron", "Theta vs P (Positron)", 100, 0, 10, 100, 0, 50);
histMPhiProton = new TH2D("MvsPhiProton", "M vs PhiProton", 120, -190, 190, 120, 0, 3);
histMPhiElectron = new TH2D("MvsPhiElectron", "M vs PhiElec", 120, -190, 190, 120, 0, 3);
histMPhiPositron = new TH2D("MvsPhiPositron", "M vs PhiPosi", 120, -190, 190, 120, 0, 3);
histMThetaElectron = new TH2D("MvsThetaElectrogdgdgn", "M vs ThetaElec", 100, 0, 50, 120, 0, 3);
histAngleMass = new TH2D("MvsThetaElectron", "M vs ThetaElec", 100, 0, 0, 120, 0, 3);
histMThetaPositron = new TH2D("MvsThetaPositron", "M vs ThetaPosi", 100, 0, 50, 120, 0, 3);
histMThetaProton = new TH2D("MvsThetaProton", "M vs ThetaProt", 100, 10, 80, 120, 0, 3);
histMVElectron = new TH2D("MvsVElectron", "M vs Vx diff pair", 100, -10, 10, 120, 0, 3);
histMVPositron = new TH2D("MvsVPositron", "M vs Vy diff pair", 100, -10, 10, 120, 0, 3);
histMV = new TH2D("MvsVProton", "M vs Vz diff pair", 100, -10, 10, 120, 0, 3);
histMCheElectron = new TH2D("MvsCheElectron", "M vs CheElec", 100, 0, 40, 120, 0, 3);
histMChePositron = new TH2D("MvsChePositron", "M vs ChePosi", 100, 0, 40, 120, 0, 3);
histMCheProton = new TH2D("MvsCheProton", "M vs CheProt", 100, 0, 40, 120, 0, 3);
histCheCheElectron = new TH2D("MvsCheElectron", "ChePosi vs CheElec", 100, 0, 40, 100, 0, 40);
histCheCheElectronInf = new TH2D("MvsCheElectrondgdgdd", "ChePosi vs CheElec", 100, 0, 40, 100, 0, 40);
histCheCheElectronSup = new TH2D("MvsCheElectronaaaa", "ChePosi vs CheElec", 100, 0, 40, 100, 0, 40);
PhiVSThetaElectron = new TH2D("PhiVSThetaelectron", "Theta vs Phi (Electron)", 100, -180, 180, 100, -0.1, 50);
PhiVSThetaPositron = new TH2D("PhiVSThetaPositron", "Theta vs Phi (Positron)", 100, -180, 180, 100, -0.1, 50);
PhiVSThetaProton = new TH2D("ThetaVSThetaProton", "Theta vs Phi (Proton)", 100, -180, 180, 100, -0.1, 70);
ThetaVSPProton = new TH2D("ThetaVSPProton", "Theta vs P (Proton)", 100, 0, 5, 100, 0, 70);
EhistMMass = new TH1F("EMMass", "Missing mass beam", 100, -1, 1);
EhistQ2 = new TH1F("EQ2", "Virtuality incoming photon (Q2)", 80, -0.1, 0.3);
EhistQ21 = new TH1F("EQ21", "Virtuality incoming photon (Q2)", 80, -0.1, 0.3);
EhistQP2 = new TH1F("EM", "Invariant Mass of the lepton pair", 120, -1, 10);
EhistM = new TH1F("EM1", "", 200, 0.0, 4.0);
EhistMAfterFid = new TH1F("EM1AfterFid", "Invariant Mass of the lepton pair AfterFid", 200, 0, 4);
EhistMprim = new TH1F("EM1prim", "Invariant Mass of the lepton pair", 200, 0, 3);
EhistEssai = new TH1F("EM2", "Invariant Mass of the lepton pair with scattered electron", 120, 0, 3);
EhistMSectorDiff = new TH2D("EM1SD", "Sector Diff vs Invariant Mass of the lepton pair", 120, 0, 3, 6, -1, 5);
EhistME = new TH2D("EMvsE", "M vs Egamma", 100, 2, 11, 100, 0, 2.5);
EhistMPElec = new TH2D("EMvsPElec", "M vs P Electron", 100, 0, 11, 100, 0, 2.5);
EhistMPPosi = new TH2D("EMvsPPosi", "M vs P Positron", 100, 0, 11, 100, 0, 2.5);
EhistMME = new TH2D("EMMvsE", "MM vs Egamma", 100, 2, 11, 100, -5, 10);
EhistMPProton = new TH2D("EMvsPProton", "M vs P Proton", 100, 0, 2, 100, 0, 2.5);
EhistMXProton = new TH2D("EMvsXProton", "M vs XProton", 100, -10, 10, 100, 0, 2.5);
EhistMRun = new TH2D("EMvsRun", "M vs RUN", 65, 4013, 4078, 30, 0, 2.5);
EhistMPhiProton = new TH2D("EMvsPhiProton", "M vs PhiProton", 120, -190, 190, 120, 0, 3);
EhistPMiss = new TH2D("EPMiss", "Direction of the missing particle", 50, -0.5, 0.5, 50, -0.5, 0.5);
EhistT = new TH1F("ET rec", "t REC", 50, -0.2, 2);
EhistMMPt = new TH2D("EMMPt", "Pt Vs MM", 100, -5, 5, 50, 0, 1);
EEnergyConservation = new TH1F("ECons. Energy", "Cons. Energy", 100, -1, 1);
EEnergyConservationEvents = new TH1F("ECons. EnergyEvents", "Cons. Energy", 100, -1, 1);
EMMQ2 = new TH2D("EMMQ2", "MM Q2", 50, -15, 15, 50, -6, 6);
EPhiVSPElectron = new TH2D("EPhiVSPelectron", "P vs Phi (Electron)", 100, -180, 180, 100, 0, 10);
EPhiVSPPositron = new TH2D("EPhiVSPPositron", "P vs Phi (Positron)", 100, -180, 180, 100, 0, 10);
EPhiVSPProton = new TH2D("EPhiVSPProton", "P vs Phi (Proton)", 100, -180, 180, 100, 0, 5);
EThetaVSPElectron = new TH2D("EThetaVSPelectron", "Theta vs P (Electron)", 100, 0, 10, 100, 0, 50);
EThetaVSPPositron = new TH2D("EThetaVSPPositron", "Theta vs P (Positron)", 100, 0, 10, 100, 0, 50);
EPhiVSThetaElectron = new TH2D("EPhiVSThetaelectron", "Theta vs Phi (Electron)", 100, -180, 180, 100, -0.1, 50);
EPhiVSThetaPositron = new TH2D("EPhiVSThetaPositron", "Theta vs Phi (Positron)", 100, -180, 180, 100, -0.1, 50);
EPhiVSThetaProton = new TH2D("EThetaVSThetaProton", "Theta vs Phi (Proton)", 100, -180, 180, 100, -0.1, 70);
EThetaVSPProton = new TH2D("EThetaVSPProton", "Theta vs P (Proton)", 100, 0, 5, 50, 0, 70);

int binPhi = 20;
int binTheta = 20;
double mintheta = 0;
double maxtheta = 140;
double minphi = -180;
double maxphi = 180;

PhiVSThetaCM = new TH2D("PhiVSThetaCMm", "Theta vs Phi (CM)", binPhi, minphi, maxphi, binTheta, mintheta, maxtheta);
PhiVSThetaCM1 = new TH2D("PhiVSThetaCM1m", "Theta vs Phi (CM)", binPhi, minphi, maxphi, binTheta, mintheta, maxtheta);
PhiVSThetaCM2 = new TH2D("PhiVSThetaCM2m", "Theta vs Phi (CM)", binPhi, minphi, maxphi, binTheta, mintheta, maxtheta);
PhiVSThetaCM3 = new TH2D("PhiVSThetaCM3m", "Theta vs Phi (CM)", binPhi, minphi, maxphi, binTheta, mintheta, maxtheta);
vertexPair = new TH1F("vertexPair", "Vertex diff pair", 50, -10, 10);
vertexElecP = new TH1F("vertexElecP", "Vertex Elec p", 50, -10, 10);
vertexPosiP = new TH1F("vertexPosP", "Vertex Posi p", 50, -10, 10);
vertexElec = new TH1F("vertexElec", "Vertex elec", 50, -10, 10);
vertexPosi = new TH1F("vertexPosi", "Vertex posi ", 50, -10, 10);
vertexProt = new TH1F("vertexProt", "Vertex proton", 50, -10, 10);
vertexProtTheta = new TH2F("vertexProtTheta", "Vertex proton vs Theta", 100, 0, 70, 50, -10, 10);
vertexProtThetaCD = new TH2F("vertexProtThetaCD", "Vertex proton vs Theta CD", 100, 0, 70, 50, -10, 10);
vertexProtThetaFD = new TH2F("vertexProtThetaFD", "Vertex proton vs Theta FD", 100, 0, 70, 50, -10, 10);
EvertexPair = new TH1F("EvertexPair", "Vertex diff pair", 50, -10, 10);
EvertexElecP = new TH1F("EvertexElecP", "Vertex Elec p", 50, -10, 10);
EvertexPosiP = new TH1F("EvertexPosP", "Vertex Posi p", 50, -10, 10);
vertexTimePair = new TH1F("vertexTimePair", "vertex time pair", 150, -2, 2);
vertexTimeP = new TH1F("vertexTimeP", "vertex time proton", 150, -2, 2);
vertexTimePP = new TH2D("vertexTimePP", "vertex time proton vs p proton", 100, 0, 3, 100, -2, 2);
Pproton = new TH1F("Pproton", "Pproton", 100, 0., 100);
EvertexTimePair = new TH1F("EvertexTimePair", "Vertex time pair", 150, -2, 2);
vertexTimePairMass = new TH2F("vertexTimePairMass", "", 100, 0, 3, 150, -2, 2);
EvertexTimeP = new TH1F("EvertexTimeP", "vertex time proton", 100, -0.5, 0.5);
EvertexTimePPCD = new TH2D("EvertexTimePP", "vertex time proton vs p proton CD", 60, 0.25, 1, 100, -0.5, 0.5);
EvertexTimePPFD = new TH2D("EvertexTimePPFD", "vertex time proton vs p proton FD", 60, 0.25, 1, 100, -0.5, 0.5);

// TCS distribution
TCSEThetaVSPElectron = new TH2D("EThetaVSPelectron", "Theta vs P (Electron)", 100, 0, 10, 100, 0, 50);
TCSEThetaVSPPositron = new TH2D("EThetaVSPPositron", "Theta vs P (Positron)", 100, 0, 10, 100, 0, 50);
TCSEThetaVSPProton = new TH2D("EThetaVSPProton", "Theta vs P (Proton)", 50, 0, 1.5, 50, 0, 70);
PPosiPElec = new TH2D("PPosiPElec", "All events;P (e+);P (e-)", 50, 0, 10, 50, 0, 10);
PPosiPElec1 = new TH2D("PPosiPElec1", "M>1.5 GeV;P (e+);P (e-)", 50, 0, 10, 50, 0, 10);

// PID histogram
SFelectron = new TH2D("SFelectron", "SF vs P (Electron)", 80, 0, 8, 100, 0., 0.35);
SFpositron = new TH2D("SFpositron", "SF vs P (Positron)", 80, 0, 8, 100, 0., 0.35);
corrSFelectron = new TH2D("corrSFelectron", "corrSF vs P (Electron)", 80, 0, 8, 100, 0.1, 0.35);
corrSFpositron = new TH2D("corrSFpositron", "corrSF vs P (Positron)", 80, 0, 8, 100, 0.1, 0.35);
ECelectron = new TH2D("ECelectron", "ECout vs ECin (Electron)", 50, 0, 0.7, 50, 0, 0.7);
ECpositron = new TH2D("ECpositron", "ECout vs ECin (Positron)", 50, 0, 0.7, 50, 0, 0.7);
EECelectron = new TH2D("EECelectron", "ECout vs ECin (Electron)", 50, 0, 0.7, 50, 0, 0.7);
EECpositron = new TH2D("EECpositron", "ECout vs ECin (Positron)", 50, 0, 0.7, 50, 0, 0.7);
CheElectron = new TH1F("Cherenkov Electron", "Cherenkov Electron", 40, 0, 40);
ChePositron = new TH1F("Cherenkov Positron", "Cherenkov Positron", 40, 0, 40);
ECheElectron = new TH1F("ECherenkov Electron", "Cherenkov Electron", 40, 0, 40);
EChePositron = new TH1F("ECherenkov Positron", "Cherenkov Positron", 40, 0, 40);
BetaProton = new TH2D("Beta P Proton", " Beta P Proton ", 100, 0, 2.5, 50, 0, 1.2);
EBetaProton = new TH2D("EBeta P Proton", " Beta P Proton ", 100, 0, 2.5, 50, 0, 1.2);

SFUelectron = new TH2D(" SF Electron vs U", " SF Electron vs U", 150, 0, 450, 50, 0., 0.35);
SFVelectron = new TH2D(" SF Electron vs V", " SF Electron vs V", 150, 0, 450, 50, 0., 0.35);
SFWelectron = new TH2D(" SF Electron vs W", " SF Electron vs W", 150, 0, 450, 50, 0., 0.35);
SFUpositron = new TH2D(" SF Positron vs U", " SF Positron vs U", 150, 0, 450, 50, 0., 0.35);
SFVpositron = new TH2D(" SF Positron vs V", " SF Positron vs V", 150, 0, 450, 50, 0., 0.35);
SFWpositron = new TH2D(" SF Positron vs W", " SF Positron vs W", 150, 0, 450, 50, 0., 0.35);

FiducialCutPCAL = new TH2D(" FiducialPCAL", "FiducialPCAL", 500, -500, 500, 500, -500, 500);
EFiducialCutPCAL = new TH2D(" EFiducialPCAL", "FiducialPCAL", 500, -500, 500, 500, -500, 500);

TH1D *VertexOtherPart = new TH1D("VertexOtherPart", "VertexOtherPart", 500, -30, 30);
TH1D *Chi2OtherPart = new TH1D("Chi2OtherPart", "Chi2OtherPart", 500, -40, 40);

EgVSPeletron = new TH2F("EgVSPeletron", "EgVSPeletron", 100, 0, 10, 100, 3, 10);
tVSPeletron = new TH2F("tVSPeletron", "tVSPeletron", 100, 0, 10, 100, 0.0, 1.);

SFMassElec = new TH2D("SFMass", "SFMassElec", 100, 0, 3, 100, 0., 0.3);
SFMassPosi = new TH2D("SFMass", "SFMassPosi", 100, 0, 3, 100, 0., 0.3);

tCheck = new TH1F("tCheck", "", 100, -1, 1);

tCheckEvents = new TH1F("tCheckEvents", "", 100, -1, 1);

tCheckBefore = new TH1F("tCheckBefore", "", 100, -1, 1);

// Check Posi PCAl/Ecal sf
SFPCALECALPosiavant = new TH2D("SFPCALECALPosi", "SFPCALECALPosi", 100, 0, 0.3, 100, 0., 0.3);
SFPCALECALPosiapres = new TH2D("SFPCALECALPosia", "SFPCALECALPosia", 100, 0, 0.3, 100, 0., 0.3);

SFPCALECALElecavant = new TH2D("SFPCALECALElec", "SFPCALECALElec", 100, 0, 0.3, 100, 0., 0.3);
SFPCALECALElecapres = new TH2D("SFPCALECALEleca", "SFPCALECALEleca", 100, 0, 0.3, 100, 0., 0.3);

SFPpcalzeroPosi = new TH2D("SFPpcalzeroPosi", "SFPpcalzeroPosi", 100, 0, 11, 100, 0., 0.3);
SFPpcalzeroElec = new TH2D("SFPpcalzeroElec", "SFPpcalzeroElec", 100, 0, 11, 100, 0., 0.3);

Q2t = new TH2F("tmass", ";t;mass", 100, 0.15, 0.8, 100, 2, 9.);

// Missing particle

PhiMissingParticle = new TH2F("PhiMissingParticle", ";phi;Pt/P", 100, -180., 180, 100, 0, 0.7);

// chi2 Part

Chi2ElectronAvant = new TH1F("Chi2ElectronAvant", "Electron (Before exclu.); #chi^{2}", 100, -6, 6);
Chi2PositronAvant = new TH1F("Chi2PositronAvant", "Positron (Before exclu.); #chi^{2}", 100, -6, 6);
Chi2ProtonAvant = new TH1F("Chi2ProtonAvant", "Proton CD (Before exclu.); #chi^{2}", 100, -10, 10);
Chi2ProtonFDAvant = new TH1F("Chi2ProtonFDAvant", "Proton FD (Before exclu.); #chi^{2}", 100, -10, 10);

Chi2Electron = new TH1F("Chi2Electron", "Electron; #chi^{2}", 100, -6, 6);
Chi2Positron = new TH1F("Chi2Positron", "Positron; #chi^{2}", 100, -6, 6);
Chi2Proton = new TH1F("Chi2Proton", "Proton CD; #chi^{2}", 100, -10, 10);
Chi2ProtonFD = new TH1F("Chi2ProtonFD", "Proton FD; #chi^{2}", 100, -10, 10);

// Plot binning
EgVST = new TH2F("EgVST", "EgVST", 80, 3, 11, 80, 0.05, 1.1);
MVST = new TH2F("MVST", "MVST", 80, 1.4, 3.1, 80, 0.05, 1.1);
MVSEg = new TH2F("MVST", "MVST", 80, 1.4, 3.1, 80, 2, 11);

Positronpt = new TH1F("PositronPt", "PositronPt", 100, 0., 4.);

PolarizationTransfer = new TH2D("PolarizationTransfer", "PolarizationTransfer", 100, 0, 11, 100, 0., 1.);
TH2D *ChePhi = new TH2D("ChePhi", "ChePhi", 200, -180, 180, 200, 0, 40);
MLPplot = new TH1F("MLPplot", "", 100, -0.05, 1.05);