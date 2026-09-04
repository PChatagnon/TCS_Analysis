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
#include "TVector3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "bib/muCLAS12class.h"
#include "bib/PhiPiEvent.h"
#include "bib/InputParser.h"

#include "hipo4/reader.h"
#include "rcdb_reader.h"

// QADB header and namespace
#include "QADB.h"
using namespace QA;

#include <ctime> // time_t
#include <cstdio>
using namespace std;

#define ADDVAR(x, name, t, tree) tree->Branch(name, x, TString(name) + TString(t))


int analysis_PhiPi()
{

	time_t begin, intermediate, end; // time_t is a datatype to store time values.

	time(&begin); // note time before execution

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
	Input input(argc, argv);

	TRandom3 *randoms = new TRandom3();

	/////////Parse command line/////////////
	bool option = input.cmdOptionExists("-option");
	isElSpectro = input.cmdOptionExists("-elSpectro");
	isGrape = input.cmdOptionExists("-grape");
	isCoincidence = input.cmdOptionExists("-coincidence");
	isCoincidence_Quasi = input.cmdOptionExists("-coincidence_quasi");
	IsInelastic = input.cmdOptionExists("-inelastic");

	cout<<"isElSpectro "<<isElSpectro<<endl;
	cout<<"isGrape "<<isGrape<<endl;
	cout<<"isCoincidence "<<isCoincidence<<endl;
	cout<<"isCoincidence_Quasi "<<isCoincidence_Quasi<<endl;
	cout<<"IsInelastic "<<IsInelastic<<endl;

	/////////////////////////////////////////

	if (input.cmdOptionExists("-energy"))
	{
		ebeam = std::stof(input.getCmdOption("-energy"));
	}

	cout<<"ebeam "<<ebeam<<endl;
	/////////End parse command line/////////////


	double nbrecEvent = 0;

	TString nameFiles = "";

	TString type = "REC";


	///////////////////////////////////////////
	// Setup the TTree output
	TString output_file = (TString)(input.getCmdOption("-o")); // argv[4]);
	TFile *outFile = new TFile(Form("output_PhiPi_%s.root", output_file.Data()), "recreate");
	
	TTree *outT = new TTree("tree", "tree");
	TTree *outT_Gen = new TTree("tree_Gen", "tree_Gen");

	TLorentzVector tree_Electron, tree_pi_plus, tree_k_plus, tree_k_minus, tree_proton, tree_Missing;
	outT->Branch("Electron", "TLorentzVector", &tree_Electron);
	outT->Branch("pi_plus", "TLorentzVector", &tree_pi_plus);
	outT->Branch("k_plus", "TLorentzVector", &tree_k_plus);
	outT->Branch("k_minus", "TLorentzVector", &tree_k_minus);
	outT->Branch("Missing", "TLorentzVector", &tree_Missing);

	int trigger_bit;
	outT->Branch("trigger_bit", &trigger_bit, "trigger_bit/I");

	std::vector<TString> fvars = {
		"evt_num",
		"run",
		"analysis_stage",
		"topology",
		
		"Q2",
		"MM",
		"M_KK",
		"M_Kppi",
		"M_Kmpi",
		"W",


	};

	std::map<TString, Float_t> outVars;
	for (size_t i = 0; i < fvars.size(); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}
	
	///////////////////////////////////////////
	


	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	int nbf = 0;
	int nbEvent = 0;
	int nbSelected = 0;
	for (Int_t i = input.getCmdIndex("-f") + 2; i < input.getCmdIndex("-ef") + 1; i++)
	{
		nbf++;
		nameFiles = TString(argv[i]);

		////////////////////////////////////////////
		// hipo reader
		hipo::reader reader;
		hipo::dictionary factory;
		hipo::event hipo_event;
		////////////////////////////////////////////
		
		reader.open(nameFiles);
		reader.readDictionary(factory);
		
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
				cout << nbEvent << " events processed in " << intermediate_time << "s" << "\n";
			}

			Y_Event ev;

			int run = 0;
			int event_nb = 0;
			
			// Get banks
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
            
			int np_input = PART.getRows();
			ev.Set_nb_part(np_input);
			
			
			
			///////////////////////////////////////////
			// Get Particles and cut on event topology
			///////////////////////////////////////////
			ev.Set_Particles(PART);

			if (!ev.pass_topology_cut())
			{
				continue;
			}
			nbSelected++;
			///////////////////////////////////////////
			// Associate detector responses and do EC cuts
			///////////////////////////////////////////
			//need to adapt this function to this final state
			//ev.Associate_detector_resp(CHE, SCIN, CALO);
			///////////////////////////////////////////

			///////////////////////////////////////////
			// Compute kinematics
			///////////////////////////////////////////
			ev.Get_Kinematics();
			///////////////////////////////////////////

			///////////////////////////////////////////
			// Fill tree
			///////////////////////////////////////////
			outVars["evt_num"] = nbEvent;
			outVars["run"] = ev.run;
			
			outVars["Q2"] = ev.Q2;
			outVars["MM"] = ev.MM;
			outVars["M_KK"] = ev.M_KK;
			outVars["M_Kppi"] = ev.M_Kppi;
			outVars["M_Kmpi"] = ev.M_Kmpi;
			outVars["W"] = ev.W;


			tree_Electron = ev.Electron.Vector;
			tree_pi_plus = ev.pi_plus.Vector;
			tree_k_minus = ev.k_minus.Vector;
			tree_k_plus = ev.k_plus.Vector;
			tree_Missing = ev.vMissing;
			
			outT->Fill();
			
		}
	}

	
	

	outFile->cd();
	outT->Write();
	outFile->Write();
	outFile->Close();

	cout << "Tree written" << endl;
	cout << "nb of file " << nbf << "\n";
	cout << "nb of events " << nbEvent << "\n";
	cout << "nb of selected " << nbSelected << "\n";

	// gROOT->ProcessLine(".q");

	time(&end); // note time after execution

	double difference = difftime(end, begin);
	printf("All this work done in only %.2lf seconds. Congratulations !\n", difference);

	gApplication->Terminate();

	return 0;
}
