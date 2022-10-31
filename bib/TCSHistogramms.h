#ifndef TCSHistogramms
#define TCSHistogramms

void init_histo();

TH1F *FinalEventst;
TH1F *FinalEventsEg;
TH1F *FinalEventsM;
TH1F *FinalEventsPhi;
TH1F *FinalEventsP;
TH1F *FinalEventsTheta;
TH2F *xihist;
TH1F *xihist1;
TH2F *thist;
TH2F *thist1;
TH2F *BeforeCuts;
TH2F *BeforeCutsPositron;
TH1F *histRECPid;
TH1F *histMMass;
TH1F *histQ2;
TH1F *histQP2;
TH1F *histM;
TH2D *histME;
TH2D *histMPElec;
TH2D *histMPPosi;
TH2D *histMME;
TH2D *histPMiss;
TH1F *histT;
TH2D *histMMPt;
TH1F *EnergyConservation;
TH2D *MMQ2;
TH2D *histQ2t;
TH2D *PhiVSPElectron;
TH2D *PhiVSPPositron;
TH2D *PhiVSPProton;
TH2D *ThetaVSPElectron;
TH2D *ThetaVSPPositron;
TH2D *histMPhiProton;
TH2D *histMPhiElectron;
TH2D *histMPhiPositron;
TH2D *histMThetaElectron;
TH2D *histAngleMass;
TH2D *histMThetaPositron;
TH2D *histMThetaProton;
TH2D *histMVElectron;
TH2D *histMVPositron;
TH2D *histMV;
TH2D *histMCheElectron;
TH2D *histMChePositron;
TH2D *histMCheProton;
TH2D *histCheCheElectron;
TH2D *histCheCheElectronInf;
TH2D *histCheCheElectronSup;
TH2D *PhiVSThetaElectron;
TH2D *PhiVSThetaPositron;
TH2D *PhiVSThetaProton;
TH2D *ThetaVSPProton;
TH1F *EhistMMass;
TH1F *EhistQ2;
TH1F *EhistQ21;
TH1F *EhistQP2;
TH1F *EhistM;
TH1F *EhistMAfterFid;
TH1F *EhistMprim;
TH1F *EhistEssai;
TH2D *EhistMSectorDiff;
TH2D *EhistME;
TH2D *EhistMPElec;
TH2D *EhistMPPosi;
TH2D *EhistMME;
TH2D *EhistMPProton;
TH2D *EhistMXProton;
TH2D *EhistMRun;
TH2D *EhistMPhiProton;
TH2D *EhistPMiss;
TH1F *EhistT;
TH2D *EhistMMPt;
TH1F *EEnergyConservation;
TH1F *EEnergyConservationEvents;
TH2D *EMMQ2;
TH2D *EPhiVSPElectron;
TH2D *EPhiVSPPositron;
TH2D *EPhiVSPProton;
TH2D *EThetaVSPElectron;
TH2D *EThetaVSPPositron;
TH2D *EPhiVSThetaElectron;
TH2D *EPhiVSThetaPositron;
TH2D *EPhiVSThetaProton;
TH2D *EThetaVSPProton;

TH2D *PhiVSThetaCM;
TH2D *PhiVSThetaCM1;
TH2D *PhiVSThetaCM2;
TH2D *PhiVSThetaCM3;
TH1F *vertexPair;
TH1F *vertexElecP;
TH1F *vertexPosiP;
TH1F *vertexElec;
TH1F *vertexPosi;
TH1F *vertexProt;
TH2F *vertexProtTheta;
TH2F *vertexProtThetaCD;
TH2F *vertexProtThetaFD;
TH1F *EvertexPair;
TH1F *EvertexElecP;
TH1F *EvertexPosiP;
TH1F *vertexTimePair;
TH1F *vertexTimeP;
TH2D *vertexTimePP;
TH1F *Pproton;
TH1F *EvertexTimePair;
TH2F *vertexTimePairMass;
TH1F *EvertexTimeP;
TH2D *EvertexTimePPCD;
TH2D *EvertexTimePPFD;

// TCS distribution
TH2D *TCSEThetaVSPElectron;
TH2D *TCSEThetaVSPPositron;
TH2D *TCSEThetaVSPProton;
TH2D *PPosiPElec;
TH2D *PPosiPElec1;

// PID histogram
TH2D *SFelectron;
TH2D *SFpositron;
TH2D *corrSFelectron;
TH2D *corrSFpositron;
TH2D *ECelectron;
TH2D *ECpositron;
TH2D *EECelectron;
TH2D *EECpositron;
TH1F *CheElectron;
TH1F *ChePositron;
TH1F *ECheElectron;
TH1F *EChePositron;
TH2D *BetaProton;
TH2D *EBetaProton;

TH2D *SFUelectron;
TH2D *SFVelectron;
TH2D *SFWelectron;
TH2D *SFUpositron;
TH2D *SFVpositron;
TH2D *SFWpositron;

TH2D *FiducialCutPCAL;
TH2D *EFiducialCutPCAL;

TH1D *VertexOtherPart;
TH1D *Chi2OtherPart;

TH2F *EgVSPeletron;
TH2F *tVSPeletron;

TH2D *SFMassElec;
TH2D *SFMassPosi;

TH1F *tCheck;

TH1F *tCheckEvents;

TH1F *tCheckBefore;

// Check Posi PCAl/Ecal sf
TH2D *SFPCALECALPosiavant;
TH2D *SFPCALECALPosiapres;

TH2D *SFPCALECALElecavant;
TH2D *SFPCALECALElecapres;

TH2D *SFPpcalzeroPosi;
TH2D *SFPpcalzeroElec;

TH2F *Q2t;

// Missing particle

TH2F *PhiMissingParticle;

// chi2 Part

TH1F *Chi2ElectronAvant;
TH1F *Chi2PositronAvant;
TH1F *Chi2ProtonAvant;
TH1F *Chi2ProtonFDAvant;

TH1F *Chi2Electron;
TH1F *Chi2Positron;
TH1F *Chi2Proton;
TH1F *Chi2ProtonFD;

// Plot binning
TH2F *EgVST;
TH2F *MVST;
TH2F *MVSEg;

TH1F *Positronpt;

TH2D *PolarizationTransfer;
TH2D *ChePhi;
TH1F *MLPplot;
#endif