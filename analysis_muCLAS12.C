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
#include "bib/muMCEvent.h"
#include "bib/muEvent.h"
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

ThetaPhi funcCM(TLorentzVector vElectron, TLorentzVector vPositron, TLorentzVector vProton)
{

	TLorentzVector vRestProton;
	TLorentzVector vPhoton;
	TLorentzVector vProtonLF;
	TLorentzVector vElectronLF;
	TLorentzVector vPositronLF;
	TLorentzVector vRestProtonLF;
	TLorentzVector vPhotonLF;
	TLorentzVector vLepton;
	TLorentzVector vIncoming;

	ThetaPhi CM;

	double phiCM;
	double thetaCM;
	double Pi = 3.14159265359;

	TLorentzVector vBeam;
	vBeam.SetPxPyPzE(0, 0, ebeam, ebeam);
	vRestProton.SetPxPyPzE(0, 0, 0, 0.938);
	vPhoton = vProton + vPositron + vElectron - vRestProton;
	phiCM = 0.0;
	thetaCM = 0.0;

	vIncoming = vProton + vPositron + vElectron;
	vPhoton.Boost(-vIncoming.BoostVector());
	vRestProton.Boost(-vIncoming.BoostVector());
	vProton.Boost(-vIncoming.BoostVector());
	vElectron.Boost(-vIncoming.BoostVector());
	vPositron.Boost(-vIncoming.BoostVector());
	if (vProton.X() == 0 && vProton.Y() == 0 && vProton.Z() == 0)
		cout << "problem1" << endl;

	vPhoton.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vElectron.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vPositron.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vRestProton.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vProton.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	// cout<<"problem proton1"<<endl;
	vElectron.RotateZ(-vRestProton.Phi());
	vPositron.RotateZ(-vRestProton.Phi());
	vProton.RotateZ(-vRestProton.Phi());
	vPhoton.RotateZ(-vRestProton.Phi());
	vRestProton.RotateZ(-vRestProton.Phi());

	phiCM = (vElectron.Phi());

	vLepton = vPositron + vElectron;
	vElectron.Boost(-vLepton.BoostVector());
	vProton.Boost(-vLepton.BoostVector());
	vPhoton.Boost(-vLepton.BoostVector());
	vRestProton.Boost(-vLepton.BoostVector());
	vPositron.Boost(-vLepton.BoostVector());

	thetaCM = vElectron.Angle(-vProton.Vect().Unit());

	CM.theta = thetaCM * TMath::RadToDeg();
	CM.phi = phiCM * TMath::RadToDeg();
	return CM;
}

int analysis_muCLAS12()
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
	TFile *outFile = new TFile(Form("output_muCLAS12_%s.root", output_file.Data()), "recreate");
	
	TTree *outT = new TTree("tree", "tree");
	TTree *outT_Gen = new TTree("tree_Gen", "tree_Gen");

	TLorentzVector tree_Electron, tree_mu_plus, tree_mu_minus, tree_Missing;
	outT->Branch("Electron", "TLorentzVector", &tree_Electron);
	outT->Branch("mu_plus", "TLorentzVector", &tree_mu_plus);
	outT->Branch("mu_minus", "TLorentzVector", &tree_mu_minus);
	outT->Branch("Missing", "TLorentzVector", &tree_Missing);

	int trigger_bit;
	outT->Branch("trigger_bit", &trigger_bit, "trigger_bit/I");

	std::vector<TString> fvars = {
		"evt_num",
		"run",
		"analysis_stage",
		"topology",
		
		"MMassProt",
		"Epho",
		"M",
		"Q2",
		"t",
		"W",

		"phi",
		"theta",

		"mu_plus_SF",
		"mu_minus_SF",	
		"n_strip_PCAL_mu_plus",
		"n_strip_PCAL_mu_minus",
		"n_strip_ECIN_mu_plus",
		"n_strip_ECIN_mu_minus",
		"n_strip_ECOUT_mu_plus",	
		"n_strip_ECOUT_mu_minus",
		"vz_elec",
		"vz_mu_plus",
		"vz_mu_minus",
	};

	std::map<TString, Float_t> outVars;
	for (size_t i = 0; i < fvars.size(); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}

	TString fvars_Gen[] = {
		"M_Gen", 
		"Q2_Gen",
		"t_Gen",
		};

	std::map<TString, Float_t> outVars_Gen;
	for (size_t i = 0; i < sizeof(fvars_Gen) / sizeof(TString); i++)
	{
		outVars_Gen[fvars_Gen[i]] = 0.;
		ADDVAR(&(outVars_Gen[fvars_Gen[i]]), fvars_Gen[i], "/F", outT_Gen);
	}

	TLorentzVector gen_Electron, gen_mu_plus, gen_mu_minus, gen_Proton;

	
	cout<<"Include all gen particles"<<endl;
	outT_Gen->Branch("gen_Electron", "TLorentzVector", &gen_Electron);
	outT_Gen->Branch("gen_mu_plus", "TLorentzVector", &gen_mu_plus);
	outT_Gen->Branch("gen_mu_minus", "TLorentzVector", &gen_mu_minus);
	outT_Gen->Branch("gen_Proton", "TLorentzVector", &gen_Proton);

	outT->Branch("gen_Electron", "TLorentzVector", &gen_Electron);
	outT->Branch("gen_mu_plus", "TLorentzVector", &gen_mu_plus);
	outT->Branch("gen_mu_minus", "TLorentzVector", &gen_mu_minus);
	outT->Branch("gen_Proton", "TLorentzVector", &gen_Proton);
	
	///////////////////////////////////////////
	


	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	int nbf = 0;
	int nbEvent = 0;
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
			if (nbEvent % 30000 == 0)
			{
				time(&intermediate);
				double intermediate_time = difftime(intermediate, begin);
				cout << nbEvent << " events processed in " << intermediate_time << "s" << "\n";
			}

			muEvent ev;
			muMCEvent MC_ev;

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
			
			MC_ev.Set_MC_Particles(MCEVENT, MCPART, isElSpectro, isGrape, isCoincidence, isCoincidence_Quasi, IsInelastic);
			MC_ev.Get_Kinematics();
			

			outVars_Gen["M_Gen"] = MC_ev.M_Gen;
			outVars_Gen["Q2_Gen"] = MC_ev.Q2_Gen;
			outVars_Gen["t_Gen"] = MC_ev.t_Gen;
			

			
			gen_Electron = MC_ev.Electron;
			gen_mu_plus = MC_ev.mu_plus;
			gen_mu_minus = MC_ev.mu_minus;
			gen_Proton = MC_ev.Proton;
			
			
            outT_Gen->Fill();
			
			///////////////////////////////////////////
			// Get Particles and cut on event topology
			///////////////////////////////////////////
			ev.Set_Particles(PART);

			if (!ev.pass_topology_cut())
			{
				continue;
			}

			///////////////////////////////////////////
			// Associate detector responses and do EC cuts
			///////////////////////////////////////////
			ev.Associate_detector_resp(CHE, SCIN, CALO);
			//ev.Associate_DC_traj(TRAJ);
			//ev.Set_Nphe_HTCC();
			///////////////////////////////////////////

			///////////////////////////////////////////
			// Compute kinematics
			///////////////////////////////////////////

			//Add the electron from MC
			double Calorimeter_resolution = 0.04;
			double smearing_factor = (1.+randoms->Gaus(0, Calorimeter_resolution)/sqrt(MC_ev.Electron.E()));
			//ev.Electron.Vector = MC_ev.Electron;
			ev.Electron.Vector.SetXYZM(smearing_factor*MC_ev.Electron.Px(), smearing_factor*MC_ev.Electron.Py(), smearing_factor*MC_ev.Electron.Pz(), me);
			ev.Correct_Momentum(MC_ev);
			ev.Get_Kinematics();


			outVars["evt_num"] = nbEvent;
			outVars["run"] = ev.run;
			outVars["MMassProt"] = ev.MMass;
			outVars["M"] = ev.M;
			outVars["Q2"] = ev.Q2;
			outVars["t"] = ev.t;
			outVars["W"] = ev.W;


			ThetaPhi cm = funcCM(ev.mu_minus.Vector, ev.mu_plus.Vector, ev.vMissing);
			outVars["phi"] = cm.phi;
			outVars["theta"] = cm.theta;

			outVars["mu_plus_SF"] =  ((ev.mu_plus.Energy(ECAL, PCAL) + ev.mu_plus.Energy(ECAL, ECIN) + ev.mu_plus.Energy(ECAL, ECOUT))) / ev.mu_plus.Vector.P();
			outVars["mu_minus_SF"] = ((ev.mu_minus.Energy(ECAL, PCAL) + ev.mu_minus.Energy(ECAL, ECIN) + ev.mu_minus.Energy(ECAL, ECOUT))) / ev.mu_minus.Vector.P();
			outVars["n_strip_PCAL_mu_plus"] = ev.mu_plus.N_strip(PCAL);
			outVars["n_strip_PCAL_mu_minus"] = ev.mu_minus.N_strip(PCAL);
			outVars["n_strip_ECIN_mu_plus"] = ev.mu_plus.N_strip(ECIN);
			outVars["n_strip_ECIN_mu_minus"] = ev.mu_minus.N_strip(ECIN);
			outVars["n_strip_ECOUT_mu_plus"] = ev.mu_plus.N_strip(ECOUT);
			outVars["n_strip_ECOUT_mu_minus"] = ev.mu_minus.N_strip(ECOUT);
			outVars["vz_elec"] = MC_ev.vz_elec_Gen;
			outVars["vz_mu_plus"] = ev.mu_plus.vertex.z;
			outVars["vz_mu_minus"] = ev.mu_minus.vertex.z;

			tree_Electron = ev.Electron.Vector;
			tree_mu_minus = ev.mu_minus.Vector;
			tree_mu_plus = ev.mu_plus.Vector;
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

	// gROOT->ProcessLine(".q");

	time(&end); // note time after execution

	double difference = difftime(end, begin);
	printf("All this work done in only %.2lf seconds. Congratulations !\n", difference);

	gApplication->Terminate();

	return 0;
}
