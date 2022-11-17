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
#include "bib/TCSEEEvent.h"
#include "bib/TCSRunSelector.h"

#include "reader.h"

#include <ctime> // time_t
#include <cstdio>
using namespace std;

#define ADDVAR(x, name, t, tree) tree->Branch(name, x, TString(name) + TString(t))

int analysisTCS_EE()
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

	TLorentzVector tree_Electron_1, tree_Electron_2, tree_Proton;
	outT->Branch("Electron_1", "TLorentzVector", &tree_Electron_1);
	outT->Branch("Electron_2", "TLorentzVector", &tree_Electron_2);
	outT->Branch("Proton", "TLorentzVector", &tree_Proton);


	TString fvars[] = {
		"Pt_Frac", "MMass", "electron_1_SF", "electron_2_SF", "electron_1_Nphe", "electron_2_Nphe",  "status_electron_1",  "status_electron_2", "status_proton"};

	std::map<TString, Float_t>
		outVars;
	for (size_t i = 0; i < sizeof(fvars) / sizeof(TString); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}

	TString nameFiles = "";
	int nbEvent = 0;
	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	for (Int_t i = 3; i < (argc - 1); i++)
	{

		if (TString(argv[i]).Contains(".hipo"))
		{
			nameFiles = TString(argv[i]);
		}
		else continue;

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

		
		while (reader.next())
		{
			nbEvent++;
			if (nbEvent % 500000 == 0)
			{
				time(&intermediate);
				double intermediate_time = difftime(intermediate, begin);

				cout << nbEvent << " events processed in " << intermediate_time << "s"
					 << "\n";
			}

			EEEvent ev;

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

			if (ev.recp == 1 && ev.recem == 2  )
			{

				if((ev.Electron[0].Energy(ECAL, PCAL) + ev.Electron[0].Energy(ECAL, ECIN))/ev.Electron[0].Vector.P()<0.2)
					continue;
				if((ev.Electron[1].Energy(ECAL, PCAL) + ev.Electron[1].Energy(ECAL, ECIN))/ev.Electron[1].Vector.P()<0.2)
					continue;
				if(ev.Proton.chi2>3.0)
					continue;
				ev.Associate_detector_resp(CHE, SCIN, CALO);
				ev.Set_Nphe_HTCC();
				ev.Compute_SF();

				outVars["Pt_Frac"] = ((ev.vBeam + ev.vRestProton - ev.Electron[0].Vector - ev.Electron[1].Vector - ev.Proton.Vector).Pt())/((ev.vBeam + ev.vRestProton - ev.Electron[0].Vector - ev.Electron[1].Vector - ev.Proton.Vector).P());
				outVars["MMass"] = (ev.vBeam + ev.vRestProton - ev.Electron[0].Vector - ev.Electron[1].Vector - ev.Proton.Vector).M2();
				outVars["electron_1_SF"] = ev.electron_1_SF;
				outVars["electron_2_SF"] = ev.electron_2_SF;
				outVars["electron_1_Nphe"] = ev.electron_1_Nphe;
				outVars["electron_2_Nphe"] = ev.electron_2_Nphe;
				outVars["status_electron_1"] = ev.Electron[0].status;
				outVars["status_electron_2"] = ev.Electron[1].status;
				outVars["status_proton"] = ev.Proton.status;

				tree_Electron_1 = ev.Electron[0].Vector;
				tree_Electron_2 = ev.Electron[1].Vector;
				tree_Proton = ev.Proton.Vector;

				outT->Fill();
			}
		}
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
