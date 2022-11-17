#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TRandom3.h"
#include "bib/TCSclass.h"
#include "bib/TCSfunc.h"
#include "bib/TCSPlotClass.h"
#include "bib/TCSMomentumCorrection.h"
#include "bib/TCSParameters.h"
#include "bib/TCSAcceptanceClass.h"
#include "bib/TCSPositronClass.h"
#include "bib/TCSFiducialCuts.h"
#include "bib/TCSBinVolumeCorrection.h"
#include "bib/TCSBSAclass.h"
#include "bib/TCSAFBclass.h"
#include "bib/TCSRRatioclass.h"
#include "bib/TCSEvent.h"
#include "bib/TCSMCEvent.h"
#include "bib/TCSRunSelector.h"

#include <ctime> // time_t
#include <cstdio>
using namespace std;

int Mix_BG()
{

	time_t begin, intermediate, end; // time_t is a datatype to store time values.

	time(&begin); // note time before execution

	Int_t argc = gApplication->Argc();
	char **argv = gApplication->Argv();

	TRandom3 *rand_gen = new TRandom3(1230);

	string inFileName = argv[3];
	TString name_mu = argv[4];
	string outFileName = inFileName.substr(inFileName.find_last_of('/') + 1, inFileName.find_last_of('.') - inFileName.find_last_of('/') - 1) + "_mixed.root";

	cout << "Check config: Input MC File " << inFileName << "Input Data File " << name_mu << " output: " << outFileName << endl;

	TChain *_ch = new TChain("tree");
	_ch->Add(name_mu);
	TTree *BG_Positron = (_ch->CopyTree("M<0.2 && Pt_Frac<0.05 && MMassBeam<0.4"))->CloneTree();

	// variables to be used to read from PU pools
	TLorentzVector *bg_Positron;
	BG_Positron->SetBranchAddress("Positron", &bg_Positron);
	cout << "done " << endl;

	TFile *oldfile = new TFile(inFileName.c_str());
	TTree *old_Events = (TTree *)oldfile->Get("tree");
	int nentries = (old_Events->GetEntriesFast());

	// list of branches to update: to do
	float input_M,input_t,input_MMassBeam,input_Pt_Frac,input_Epho,input_qp2;
	TLorentzVector *in_Positron = 0;
	TLorentzVector *in_Electron = 0;
	TLorentzVector *in_Proton = 0;
	old_Events->SetBranchAddress("Positron", &in_Positron);
	old_Events->SetBranchAddress("Electron", &in_Electron);
	old_Events->SetBranchAddress("Proton", &in_Proton);
	old_Events->SetBranchAddress("M", &input_M);
	old_Events->SetBranchAddress("t", &input_t);
	old_Events->SetBranchAddress("MMassBeam", &input_MMassBeam);
	old_Events->SetBranchAddress("Pt_Frac", &input_Pt_Frac);
	old_Events->SetBranchAddress("Epho", &input_Epho);
	old_Events->SetBranchAddress("qp2", &input_qp2);
	//old_Events->SetBranchAddress("xi", &input_xi);
	//old_Events->SetBranchAddress("s", &input_s);
	//old_Events->SetBranchAddress("L", &input_L);
	//old_Events->SetBranchAddress("L0", &input_L0);
	//old_Events->SetBranchAddress("theta", &input_theta);
	//old_Events->SetBranchAddress("phi", &input_phi);
	//old_Events->SetBranchAddress("positron_SF", &input_positron_SF);
	old_Events->SetBranchStatus("*", 1); // activate all branches to copy

	// Open output mixed file
	TFile *fOutMixed = new TFile(outFileName.c_str(), "RECREATE");
	if (!fOutMixed->IsOpen())
	{
		cout << "Error opening output file. Abort." << endl;
		return 5;
	}

	TTree *mixed_Events = old_Events->CloneTree(0);

	int nbEvent = 0;

	Event ev;

	while (nbEvent < nentries /*&& nbEvent < 10*/)
	{
		old_Events->GetEntry(nbEvent);
		nbEvent++;

		int i_event = rand_gen->Rndm() * BG_Positron->GetEntries();
		BG_Positron->GetEntry(i_event);

		in_Positron=bg_Positron;

		ev.Set_Vectors(*in_Electron, *in_Positron, *in_Proton);
		ev.Get_Kinematics();

		
		input_M=ev.M;
		input_t=ev.t;
		input_MMassBeam=ev.MMassBeam;
		input_Pt_Frac=ev.Pt_Frac;
		input_Epho=ev.Epho;
		input_qp2=ev.qp2;

		mixed_Events->Fill();

		if (nbEvent % 30000 == 0)
		{
			time(&intermediate);
			double intermediate_time = difftime(intermediate, begin);

			cout << nbEvent << " events processed in " << intermediate_time << "s"
				 << "\n";
		}
	}

	// Write mixed tree to file and close everything
	fOutMixed->cd();
	cout << "Writes " << fOutMixed->GetName() << endl;
	mixed_Events->Write();
	fOutMixed->Close();

	time(&end); // note time after execution

	double difference = difftime(end, begin);
	printf("All this work done in only %.2lf seconds. Congratulations !\n", difference);

	gApplication->Terminate();

	return 0;
}
