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
#include "TH3F.h"
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
#include "bib/TCSBGEvent.h"
#include "bib/TCSRunSelector.h"

#include "reader.h"

#include <ctime> // time_t
#include <cstdio>
using namespace std;

#define ADDVAR(x, name, t, tree) tree->Branch(name, x, TString(name) + TString(t))

int analysisTCS_BG()
{

	time_t begin, intermediate, end;

	time(&begin);

	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(111);
	gStyle->SetPalette(55);
	gStyle->SetLabelSize(.05, "xyz");
	gStyle->SetTitleSize(.05, "xyz");
	gStyle->SetTitleSize(.07, "t");
	gStyle->SetMarkerStyle(13);
	gStyle->SetOptFit(1);

	Int_t argc = gApplication->Argc();
	char **argv = gApplication->Argv();

	TString output_file = (TString)(argv[argc - 1]);
	TFile *outFile = new TFile(Form("outputTCS_%s.root", output_file.Data()), "recreate");
	TTree *outT = new TTree("tree", "tree");

	TLorentzVector tree_Electron, tree_Positron, tree_Proton, tree_PionM;
	outT->Branch("Electron", "TLorentzVector", &tree_Electron);
	outT->Branch("Positron", "TLorentzVector", &tree_Positron);
	outT->Branch("Proton", "TLorentzVector", &tree_Proton);
	outT->Branch("PionM", "TLorentzVector", &tree_PionM);

	TString fvars[] = {
		"MMass", "MMassProton", "positron_SF", "electron_SF", "positron_Nphe", "electron_Nphe"};

	std::map<TString, Float_t>
		outVars;
	for (size_t i = 0; i < sizeof(fvars) / sizeof(TString); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}
	
	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	TString nameFiles = TString(argv[3]);

	////////////////////////////////////////////
	// hipo reader
	hipo::reader reader;
	hipo::dictionary factory;
	hipo::event hipo_event;
	////////////////////////////////////////////

	reader.open(nameFiles);
	reader.readDictionary(factory);
	factory.show();

	hipo::bank EVENT(factory.getSchema("REC::Event"));
	hipo::bank PART(factory.getSchema("REC::Particle"));
	hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
	hipo::bank CHE(factory.getSchema("REC::Cherenkov"));
	hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
	hipo::bank RUN(factory.getSchema("RUN::config"));
	hipo::bank MCPART(factory.getSchema("MC::Particle"));
	hipo::bank MCEVENT(factory.getSchema("MC::Event"));
	hipo::bank TRACK(factory.getSchema("REC::Track"));
	hipo::bank TRAJ(factory.getSchema("REC::Traj"));

	outFile->cd();

	int nbEvent = 0;
	while (reader.next())
	{
		nbEvent++;
		if (nbEvent % 30000 == 0)
		{
			time(&intermediate);
			double intermediate_time = difftime(intermediate, begin);

			cout << nbEvent << " events processed in " << intermediate_time << "s"
				 << "\n";
		}

		BGEvent ev;

		reader.read(hipo_event);
		hipo_event.getStructure(MCPART);
		hipo_event.getStructure(MCEVENT);
		hipo_event.getStructure(RUN);
		hipo_event.getStructure(PART);
		hipo_event.getStructure(SCIN);
		hipo_event.getStructure(CHE);
		hipo_event.getStructure(CALO);
		hipo_event.getStructure(EVENT);
		hipo_event.getStructure(TRAJ);
		hipo_event.getStructure(TRACK);

		if (PART.getSize() < 1)
			continue;

		ev.Set_Particles(PART);
		ev.Apply_EC_Cuts(CALO);
		ev.Associate_detector_resp(CHE, SCIN);
		ev.Set_Nphe_HTCC();
		ev.Compute_SF();

		outVars["MMass"] = (ev.vBeam+ev.vRestProton-ev.Positron.Vector-ev.Electron.Vector-ev.Proton.Vector-ev.PionM.Vector).M2();
		outVars["MMassProton"] = (ev.vBeam+ev.vRestProton-ev.Positron.Vector-ev.Electron.Vector-ev.PionM.Vector).M2();
		outVars["positron_SF"] = ev.positron_SF;
		outVars["electron_SF"] = ev.electron_SF;
		outVars["positron_Nphe"] = ev.positron_Nphe;
		outVars["electron_Nphe"] = ev.electron_Nphe;

		tree_Electron = ev.Electron.Vector;
		tree_Positron = ev.Positron.Vector;
		tree_Proton = ev.Proton.Vector;
		tree_PionM = ev.PionM.Vector;

		outT->Fill();
	}

	outT->Write();
	outFile->Write();
	outFile->Close();

	time(&end);

	double difference = difftime(end, begin);
	printf("All this work done in only %.2lf seconds. Congratulations !\n", difference);

	gApplication->Terminate();

	return 0;
}
