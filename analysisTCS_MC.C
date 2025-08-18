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

int analysisTCS_MC()
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

	Parameters InputParameters("InputOptions/InputOptionsTCS_Standard.txt");

	Int_t argc = gApplication->Argc();
	char **argv = gApplication->Argv();
	Input input(argc, argv);

	/////////Parse command line/////////////
	IsEE_BG = input.cmdOptionExists("-IsEE_BG");
	IsTCSGen = input.cmdOptionExists("-IsTCSGen");
	IsGrape = input.cmdOptionExists("-IsGrape");
	IsJPsi = input.cmdOptionExists("-IsJPsi");
	IsElSpectro = input.cmdOptionExists("-IsElSpectro");
	RGA_Fall2018 = input.cmdOptionExists("-RGA_Fall2018");
	RGA_Spring2019 = input.cmdOptionExists("-RGA_Spring2019");
	inbending = !input.cmdOptionExists("-outbending");
	PCAL_study = input.cmdOptionExists("-PCAL");
	CALO_study = input.cmdOptionExists("-CALO");
	Lepton_ID_check = input.cmdOptionExists("-Lepton_ID_check");
	DC_Traj_check = input.cmdOptionExists("-DC_Traj_check");
	all_Gen_vector = input.cmdOptionExists("-all_Gen_vector");
	QA_Golden = input.cmdOptionExists("-QA_Golden");
	no_QADB = input.cmdOptionExists("-no_QADB");
	no_RCDB = input.cmdOptionExists("-no_RCDB");
	inclusive_topology = input.cmdOptionExists("-inclusive_topology");
	FD_mom_corr = input.cmdOptionExists("-FD_mom_corr");
	v2_mom_corr = input.cmdOptionExists("-v2_mom_corr");
	/////////////////////////////////////////

	if (input.cmdOptionExists("-energy"))
	{
		ebeam = std::stof(input.getCmdOption("-energy"));
	}

	if (input.cmdOptionExists("-usage"))
	{
		cout << "Use as : clas12root -l analysisTCS_MC.C -a NewacceptanceTCS_newSimuLargeStats.root -o ouputname -f files -ef -inbending\n";
		cout << "Options available:\n";
		cout << "-IsEE_BG\n";
		cout << "-IsTCSGen\n";
		cout << "-IsGrape\n";
		cout << "-IsJPsi\n";
		cout << "-RGA_Fall2018\n";
		cout << "-RGA_Spring2019\n";
		cout << "-outbending\n";
		cout << "-PCAL\n";
		cout << "-Lepton_ID_check\n";
		cout << "-DC_Traj_check\n";
		cout << "-all_Gen_vector\n";
		cout << "-energy\n";

		gApplication->Terminate();
	}
	/////////End parse command line/////////////

	cout << "////////////////////////////////////////////"
		 << "\n";
	cout << "Run with the following options : "
		 << "\n";
	cout << "Inbending : " << inbending << "\n";
	cout << "RGA_Fall2018 : " << RGA_Fall2018 << "\n";
	cout << "RGA_Spring2019 : " << RGA_Spring2019 << "\n";
	cout << "IsTCSGen : " << IsTCSGen << "\n";
	cout << "IsGrape : " << IsGrape << "\n";
	cout << "IsJPsi : " << IsJPsi << "\n";
	cout << "IsElSpectro : " << IsElSpectro << "\n";
	cout << "IsEE_BG : " << IsEE_BG << "\n";
	cout << "Lepton ID check : " << Lepton_ID_check << "\n";
	cout << "PCAL var. : " << PCAL_study << "\n";
	cout << "CALO study : " << CALO_study << "\n";
	cout << "Run with energy " << ebeam << " GeV\n";
	cout << "////////////////////////////////////////////"
		 << "\n";

	/////////Instanciate QADB///////////
	QADB *qa = new QADB("pass2");
	////////////////////////////////////

	double nbrecEvent = 0;
	int nbf = 0;
	double nEventTCS = 0;
	double denom = 0;
	double nCD = 0;
	double nFD = 0;

	int nbevent_after_EC = 0;
	int nbevent_after_posi = 0;

	int AfterCuts = 0;

	TString nameFiles = "";

	TString type = "REC";

	/////////////////////////////////////////////
	// RCDB setup
	/////////////////////////////////////////////
	rcdb_root rcdb("InputOptions/rcdb.root");

	//////////////////////////////////////////////
	// Acceptance setup
	//////////////////////////////////////////////

	// Acceptance Acc_TCS(TString(argv[3]), 4, 3, 3, 36, 13);
	Acceptance Acc_TCS(TString(input.getCmdOption("-a")), 4, 3, 3, 36, 13);

	Acc_TCS.Draw_Acc();
	Acc_TCS.Draw_Error();

	////////////////////////////////////////////
	// Init run selector
	////////////////////////////////////////////

	RunSelector Run_Selector;

	////////////////////////////////////////////
	// Bin volume correction
	////////////////////////////////////////////

	BinVolumeCorrection Bin_Volume_Corr(Acc_TCS);

	//////////////////////////////////////////////
	// Momentum Correction
	//////////////////////////////////////////////
	MomentumCorrection MomCorr;
	Energy_loss EnergyLoss(inbending, RGA_Fall2018);
	Momentum_Corrections_Richard Momentum_Correction(inbending, v2_mom_corr, FD_mom_corr);

	////////////////////////////////////////////
	// Instanciate observables
	////////////////////////////////////////////

	float t1_f[5] = {0.15, 0.25, 0.34, 0.48, 0.8};
	float t1b_f[5] = {0.15, 0.35, 0.45, 0.55, 0.8};
	float Xi_f[4] = {0.0, 0.12, 0.15, 0.4};
	float M_f[5] = {1.5, 1.7, 2., 2.5, 3.};
	float Eg_f[4] = {4., 6., 8., 10.};

	TCSBSA BSA_vs_t("t", 4, t1_f, 10);
	TCSBSA BSA_vs_Xi("xi", 3, Xi_f, 10);
	TCSBSA BSA_vs_M("M", 4, M_f, 10);

	TCSAFB AFB_vs_t("t", 4, t1_f);
	TCSAFB AFB_vs_t1("t1", 4, t1_f);
	TCSAFB AFB_vs_t2("t2", 4, t1b_f);
	TCSAFB AFB_vs_M("M", 4, M_f);
	TCSAFB AFB_vs_Eg("Eg", 3, Eg_f);

	TCSRRatio RRatio_vs_t("t", 4, t1_f);
	TCSRRatio RRatio_vs_Xi("Xi", 3, Xi_f);

	///////////////////////////////////////////
	// Setup the TTree output
	TString output_file = (TString)(input.getCmdOption("-o")); // argv[4]);
	TFile *outFile = new TFile(Form("outputTCS_%s.root", output_file.Data()), "recreate");
	// TFile *outFile = new TFile("outputTCS_"+output_file+".root", "recreate");
	TTree *outT = new TTree("tree", "tree");
	TTree *outT_Gen = new TTree("tree_Gen", "tree_Gen");

	TLorentzVector tree_Electron, tree_Positron, tree_Proton, tree_Missing;
	outT->Branch("Electron", "TLorentzVector", &tree_Electron);
	outT->Branch("Positron", "TLorentzVector", &tree_Positron);
	outT->Branch("Proton", "TLorentzVector", &tree_Proton);
	outT->Branch("Missing", "TLorentzVector", &tree_Missing);

	int trigger_bit;
	outT->Branch("trigger_bit", &trigger_bit, "trigger_bit/I");

	std::vector<TString> fvars = {
		"evt_num",
		"topology_FS",
		"t",
		"t_min",
		"MMassBeam",
		"Epho",
		"Epho_Gen",
		"qp2",
		"M",
		"xi",
		"s",
		"L",
		"L0",
		"Pt_Frac",
		"Q2",
		"Q2_Gen",
		"theta",
		"phi",
		"theta_GJ_Gen",
		"phi_GJ_Gen",
		"helicity",
		"polaT",
		"positron_SF",
		"electron_SF",
		"positron_score",
		"electron_score",
		"weight",
		"acc",
		"acc_error",
		"real_flux_Gen",
		"virtual_flux_Gen",
		"virtual_flux_Frixione_Gen",
		"run",
		"analysis_stage",
		"topology",
		"positron_Nphe",
		"electron_Nphe",
		"positron_HTCCt",
		"electron_HTCCt",
		"pass_EC_cut",
		"positron_HTCC_ECAL_match",
		"electron_HTCC_ECAL_match",
		"status_elec",
		"status_posi",
		"status_prot",
		"vx_elec",
		"vy_elec",
		"vz_elec",
		"vx_posi",
		"vy_posi",
		"vz_posi",
		"vx_prot",
		"vy_prot",
		"vz_prot",
		"chi2_proton",
		"PCAL_sector_elec",
		"PCAL_sector_posi",
		"lead_lep_p",
		"sub_lead_lep_p",
		"lead_lep_theta",
		"sub_lead_lep_theta",
		"Triangular_Cut_elec",
		"Triangular_Cut_posi",
		"CM_gamma_energy",
		"CM_gamma_energy_2",
		"Q2_true",
		"E_k",
		"E_k_2",
		"PCAL_x_elec",
		"PCAL_y_elec",
		"PCAL_x_posi",
		"PCAL_y_posi",
	};

	if (PCAL_study)
	{
		fvars.insert(fvars.end(), {"PCAL_x_elec",
								   "PCAL_y_elec",
								   "PCAL_U_elec",
								   "PCAL_V_elec",
								   "PCAL_W_elec",
								   "PCAL_x_posi",
								   "PCAL_y_posi",
								   "PCAL_U_posi",
								   "PCAL_V_posi",
								   "PCAL_W_posi"});
	}

	if (CALO_study)
	{
		fvars.insert(fvars.end(), {"PCAL_x_elec_rot", "PCAL_y_elec_rot", "PCAL_z_elec_rot",
								   "PCAL_x_posi_rot", "PCAL_y_posi_rot", "PCAL_z_posi_rot",
								   "PCAL_hx_elec_rot", "PCAL_hy_elec_rot", "PCAL_hz_elec_rot",
								   "PCAL_hx_posi_rot", "PCAL_hy_posi_rot", "PCAL_hz_posi_rot",
								   "ECIN_x_elec_rot", "ECIN_y_elec_rot", "ECIN_z_elec_rot",
								   "ECIN_x_posi_rot", "ECIN_y_posi_rot", "ECIN_z_posi_rot",
								   "ECIN_hx_elec_rot", "ECIN_hy_elec_rot", "ECIN_hz_elec_rot",
								   "ECIN_hx_posi_rot", "ECIN_hy_posi_rot", "ECIN_hz_posi_rot"});
	}

	if (Lepton_ID_check)
	{
		fvars.insert(fvars.end(), {"SFPCAL_elec", "SFECIN_elec", "SFECOUT_elec",
								   "SFPCAL_posi", "SFECIN_posi", "SFECOUT_posi",
								   "M2PCAL_elec", "M2ECIN_elec", "M2ECOUT_elec",
								   "M2PCAL_posi", "M2ECIN_posi", "M2ECOUT_posi"});
	}
	if (DC_Traj_check)
	{
		fvars.insert(fvars.end(), {"DC_R1_elec_x", "DC_R1_elec_y", "DC_R1_elec_z",
								   "DC_R1_posi_x", "DC_R1_posi_y", "DC_R1_posi_z",
								   "DC_R1_prot_x", "DC_R1_prot_y", "DC_R1_prot_z",

								   "DC_R2_elec_x", "DC_R2_elec_y", "DC_R2_elec_z",
								   "DC_R2_posi_x", "DC_R2_posi_y", "DC_R2_posi_z",
								   "DC_R2_prot_x", "DC_R2_prot_y", "DC_R2_prot_z",

								   "DC_R3_elec_x", "DC_R3_elec_y", "DC_R3_elec_z",
								   "DC_R3_posi_x", "DC_R3_posi_y", "DC_R3_posi_z",
								   "DC_R3_prot_x", "DC_R3_prot_y", "DC_R3_prot_z"});
	}

	std::map<TString, Float_t> outVars;
	for (size_t i = 0; i < fvars.size(); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}

	TString fvars_Gen[] = {
		"weight", "evt_num", "t_Gen", "t_min_Gen", "MMassBeam_Gen", "Epho_Gen", "qp2_Gen", "M_Gen_1", "M_Gen_2", "Pt_Frac_Gen", "Q2_Gen",
		"vz_elec_Gen", "vz_posi_Gen", "vz_prot_Gen",
		"theta_GJ_Gen","phi_GJ_Gen",
		"theta_Gen", "phi_Gen", "real_flux_Gen", "virtual_flux_Gen", "virtual_flux_Frixione_Gen"};

	std::map<TString, Float_t> outVars_Gen;
	if (IsGrape || IsTCSGen || IsJPsi || IsElSpectro)
	{
		for (size_t i = 0; i < sizeof(fvars_Gen) / sizeof(TString); i++)
		{
			outVars_Gen[fvars_Gen[i]] = 0.;
			ADDVAR(&(outVars_Gen[fvars_Gen[i]]), fvars_Gen[i], "/F", outT_Gen);
		}
	}

	TLorentzVector gen_Electron, gen_Positron, gen_Proton;

	if (all_Gen_vector)
	{
		outT_Gen->Branch("gen_Electron", "TLorentzVector", &gen_Electron);
		outT_Gen->Branch("gen_Positron", "TLorentzVector", &gen_Positron);
		outT_Gen->Branch("gen_Proton", "TLorentzVector", &gen_Proton);

		outT->Branch("gen_Electron", "TLorentzVector", &gen_Electron);
		outT->Branch("gen_Positron", "TLorentzVector", &gen_Positron);
		outT->Branch("gen_Proton", "TLorentzVector", &gen_Proton);
	}
	///////////////////////////////////////////

	int nbtc = 0;
	int nbJPSI = 0;

	///////////////////////////////////////////
	// TMVA PID for Positron
	///////////////////////////////////////////

	// string path_ML_weights = "ML_weights";
	string path_ML_weights = "ML_weights_pass2";

	TString positron_bdt_weights;
	TString electron_bdt_weights;

	if (inbending && RGA_Fall2018) /// Inbending Fall 2018
	{

		cout << "////////////////////////////////////////////"
			 << "\n";
		cout << "Weight for Inbending Fall 2018 " << ebeam << " GeV\n";
		cout << "////////////////////////////////////////////"
			 << "\n";

		if (IsEE_BG)
		{
			positron_bdt_weights = path_ML_weights + "/F18inneg/TMVAClassification_BDT.weights.xml";
			electron_bdt_weights = path_ML_weights + "/F18inneg/TMVAClassification_BDT.weights.xml";
		}
		else
		{
			positron_bdt_weights = path_ML_weights + "/F18inpos/TMVAClassification_BDT.weights.xml";
			electron_bdt_weights = path_ML_weights + "/F18inneg/TMVAClassification_BDT.weights.xml";
		}
	}
	else if (inbending && RGA_Spring2019) /// Inbending Spring 2019
	{

		cout << "////////////////////////////////////////////"
			 << "\n";
		cout << "Weight for Inbending Spring 2019 " << ebeam << " GeV\n";
		cout << "////////////////////////////////////////////"
			 << "\n";

		if (IsEE_BG)
		{
			positron_bdt_weights = path_ML_weights + "/S19neg/TMVAClassification_BDT.weights.xml";
			electron_bdt_weights = path_ML_weights + "/S19neg/TMVAClassification_BDT.weights.xml";
		}
		else
		{
			positron_bdt_weights = path_ML_weights + "/S19pos/TMVAClassification_BDT.weights.xml";
			electron_bdt_weights = path_ML_weights + "/S19neg/TMVAClassification_BDT.weights.xml";
		}
	}
	else if (!inbending) /// Outbending  Fall 2018
	{

		cout << "////////////////////////////////////////////"
			 << "\n";
		cout << "Weight for Outbending Fall 2018 " << ebeam << " GeV\n";
		cout << "////////////////////////////////////////////"
			 << "\n";

		if (IsEE_BG)
		{
			positron_bdt_weights = path_ML_weights + "/F18outneg/TMVAClassification_BDT.weights.xml";
			electron_bdt_weights = path_ML_weights + "/F18outneg/TMVAClassification_BDT.weights.xml";
		}
		else
		{
			positron_bdt_weights = path_ML_weights + "/F18outpos/TMVAClassification_BDT.weights.xml";
			electron_bdt_weights = path_ML_weights + "/F18outneg/TMVAClassification_BDT.weights.xml";
		}
	}

	PositronIdentification PositronPID("BDT", positron_bdt_weights, 0.0, 4.0);
	PositronPID.InitializeBDT_new_PositronIdentification();

	PositronIdentification ElectronPID("BDT", electron_bdt_weights, 0.0, 4.0);
	ElectronPID.InitializeBDT_new_PositronIdentification();

	///////////////////////////////////////////
	// Plots
	///////////////////////////////////////////
	TCSPlots Plots;
	Plots.Initialize_1D();
	Plots.Initialize_2D();
	Plots.SetOutputFolder("Plots");

	// chi2cut proton
	double meanFD = 0.26;
	double meanCD = 0.81;
	double sigmaFD = 1.207;
	double sigmaCD = 1.972;

	int seed = 0;
	TRandom *polarizationGene = new TRandom(seed);
	TRandom3 *ChoiceEvent = new TRandom3(seed);

	int corrrad = 0;
	int nbEvent = 0;

	////////////////////////////////////////////
	// Initialize RCDB flags
	bool RCDB_read = false;
	double beam_current = 0.0;
	string beam_current_requested = "";
	int run_data = 0;
	////////////////////////////////////////////

	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	// for (Int_t i = 6; i < (argc); i++)
	// cout<<input.getCmdIndex("-f")<<"\n";
	// cout<<input.getCmdIndex("-ef")<<"\n";
	for (Int_t i = input.getCmdIndex("-f") + 2; i < input.getCmdIndex("-ef") + 1; i++)
	{
		if (TString(argv[i]).Contains("MC") || IsGrape || IsTCSGen || IsJPsi || IsElSpectro)
		{
			IsData = false;
		}

		if (TString(argv[i]).Contains(".hipo"))
		{
			nbf++;
			nameFiles = TString(argv[i]);
		}

		if (TString(argv[5]).Contains(".root"))
		{
			IsHipo = false;
			nameFiles = TString(argv[i]);
		}

		////////////////////////////////////////////
		cout << "////////////////////////////////////////////"
			 << "\n";
		if (IsData)
			cout << "Running on Data"
				 << "\n";
		else if (IsTCSGen)
			cout << "Running on TCSGen Simulation"
				 << "\n";
		else if (IsGrape)
			cout << "Running on Grape Simulation"
				 << "\n";
		else if (IsJPsi)
			cout << "Running on JPsi Simulation"
				 << "\n";
		else if (IsElSpectro)
			cout << "Running on ElSpectro Simulation"
				 << "\n";

		cout << TString(argv[i]) << "\n";
		cout << "Is hipo ? " << IsHipo << "\n";
		cout << "////////////////////////////////////////////"
			 << "\n";
		cout << "Run with the following options : "
			 << "\n";
		cout << "Inbending : " << inbending << "\n";
		cout << "RGA_Fall2018 : " << RGA_Fall2018 << "\n";
		cout << "IsTCSGen : " << IsTCSGen << "\n";
		cout << "IsGrape : " << IsGrape << "\n";
		cout << "IsJPsi : " << IsJPsi << "\n";
		cout << "IsElSpectro : " << IsElSpectro << "\n";
		cout << "IsEE_BG : " << IsEE_BG << "\n";
		cout << "////////////////////////////////////////////"
			 << "\n";
		////////////////////////////////////////////

		////////////////////////////////////////////
		// hipo reader
		hipo::reader reader;
		hipo::dictionary factory;
		hipo::event hipo_event;
		////////////////////////////////////////////

		////////////////////////////////////////////
		// Root TTree reader
		TFile *input_root_file;
		TTree *input_tree;
		int nentries;
		TLorentzVector *input_Electron = 0;
		TLorentzVector *input_Positron = 0;
		TLorentzVector *input_Proton = 0;
		float input_t, input_MMassBeam, input_Epho, input_qp2, input_M, input_xi, input_Pt_Frac, input_theta, input_phi, input_positron_SF, input_electron_SF, input_weight;
		float input_run, input_acc, input_acc_error, input_real_flux, input_virtual_flux;
		float input_s, input_L0, input_L;
		////////////////////////////////////////////

		if (IsHipo)
		{
			reader.open(nameFiles);
			reader.readDictionary(factory);
			// factory.show();
		}

		if (!IsHipo)
		{
			input_root_file = new TFile(nameFiles);
			input_tree = (TTree *)input_root_file->Get("tree");
			nentries = (input_tree->GetEntriesFast());
			input_tree->SetBranchAddress("Electron", &input_Electron);
			input_tree->SetBranchAddress("Positron", &input_Positron);
			input_tree->SetBranchAddress("Proton", &input_Proton);
			input_tree->SetBranchAddress("t", &input_t);
			input_tree->SetBranchAddress("MMassBeam", &input_MMassBeam);
			input_tree->SetBranchAddress("Epho", &input_Epho);
			input_tree->SetBranchAddress("qp2", &input_qp2);
			input_tree->SetBranchAddress("M", &input_M);
			input_tree->SetBranchAddress("xi", &input_xi);
			input_tree->SetBranchAddress("s", &input_s);
			input_tree->SetBranchAddress("L", &input_L);
			input_tree->SetBranchAddress("L0", &input_L0);
			input_tree->SetBranchAddress("Pt_Frac", &input_Pt_Frac);
			input_tree->SetBranchAddress("theta", &input_theta);
			input_tree->SetBranchAddress("phi", &input_phi);
			input_tree->SetBranchAddress("positron_SF", &input_positron_SF);
			input_tree->SetBranchAddress("electron_SF", &input_electron_SF);
			input_tree->SetBranchAddress("weight", &input_weight);
			input_tree->SetBranchAddress("run", &input_run);
			input_tree->SetBranchAddress("acc", &input_acc);
			input_tree->SetBranchAddress("acc_error", &input_acc_error);
			input_tree->SetBranchAddress("real_flux", &input_real_flux);
			input_tree->SetBranchAddress("virtual_flux", &input_virtual_flux);
		}

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

		while (((reader.next() && IsHipo) || (nbEvent < nentries && !IsHipo)) /*&& nbEvent < 100000*/)
		{

			nbEvent++;
			if (nbEvent % 30000 == 0)
			{
				time(&intermediate);
				double intermediate_time = difftime(intermediate, begin);

				cout << nbEvent << " events processed in " << intermediate_time << "s"
					 << "\n";
			}

			// cout << "event" << endl;
			double w = 1; // MCfluxBH[0]*MCpsfBH[0]*MCcsBH[0];
			int polarization;
			Event ev;
			MCEvent MC_ev;

			int run = 0;
			int event_nb = 0;
			// int trigger_bit = 0;

			if (IsHipo)
			{

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

				if (MCPART.getSize() < 1 && (!IsData))
					continue;

				// Number of total event
				Plots.Fill_1D("evt_count", 0, 1);

				run = RUN.getInt("run", 0);
				event_nb = RUN.getInt("event", 0);
				trigger_bit = RUN.getLong("trigger", 0);
				int np_input = PART.getRows();
				ev.Set_nb_part(np_input);
				ev.Set_trigger_bit(trigger_bit);

				if (!IsData)
				{

					MC_ev.Set_MC_Particles(MCEVENT, MCPART, IsGrape, IsJPsi, IsElSpectro, IsTCSGen);
					MC_ev.Get_Kinematics(IsGrape);

					if (IsTCSGen || IsJPsi)
						w = MC_ev.w;

					/*if (IsJPsi)
					{
						float MC_factor_1 = MCEVENT.getFloat("ptarget", 0);
						float MC_factor_2 = MCEVENT.getFloat("pbeam", 0);
						float MC_factor_3 = MCEVENT.getFloat("ebeam", 0);
						w = MC_factor_1 * MC_factor_2 * MC_factor_3;
					}*/

					ev.Set_Weight(w);

					outVars_Gen["t_Gen"] = MC_ev.t_Gen;
					outVars_Gen["t_min_Gen"] = MC_ev.t_min_Gen;
					outVars_Gen["MMassBeam_Gen"] = MC_ev.MMassBeam_Gen;
					outVars_Gen["Epho_Gen"] = MC_ev.Epho_Gen;
					outVars_Gen["qp2_Gen"] = MC_ev.qp2_Gen;
					outVars_Gen["M_Gen_1"] = MC_ev.M_Gen_1;
					outVars_Gen["M_Gen_2"] = MC_ev.M_Gen_2;
					outVars_Gen["Pt_Frac_Gen"] = MC_ev.Pt_Frac_Gen;
					outVars_Gen["Q2_Gen"] = MC_ev.Q2_Gen;
					outVars_Gen["theta_Gen"] = MC_ev.theta_Gen;
					outVars_Gen["phi_Gen"] = MC_ev.phi_Gen;
					outVars_Gen["theta_GJ_Gen"] = MC_ev.theta_GJ_Gen;
					outVars_Gen["phi_GJ_Gen"] = MC_ev.phi_GJ_Gen;
					outVars_Gen["real_flux_Gen"] = MC_ev.real_flux_Gen;
					outVars_Gen["virtual_flux_Gen"] = MC_ev.virtual_flux_Gen;
					outVars_Gen["virtual_flux_Frixione_Gen"] = MC_ev.virtual_flux_Frixione_Gen;
					outVars_Gen["evt_num"] = nbEvent;
					outVars_Gen["weight"] = w;
					outVars_Gen["vz_elec_Gen"] = MC_ev.vz_elec_Gen;
					outVars_Gen["vz_posi_Gen"] = MC_ev.vz_posi_Gen;
					outVars_Gen["vz_prot_Gen"] = MC_ev.vz_prot_Gen;

					if (all_Gen_vector)
					{
						gen_Electron = MC_ev.Electron_2;
						gen_Positron = MC_ev.Positron;
						gen_Proton = MC_ev.Proton;
					}

					outT_Gen->Fill();
				}

				///////////////////////////////////////////
				// Filter good runs for data only using QADB
				///////////////////////////////////////////
				// if (!Run_Selector.Is_Good_Run(run) && IsData && RGA_Fall2018)
				bool Keep_event = true;
				if (IsData && !no_QADB)
				{
					Keep_event = qa->OkForAsymmetry(run, event_nb);
					if (QA_Golden)
					{
						Keep_event = qa->Golden(run, event_nb);
					}
					// if (no_QADB)
					//{
					//	qa->SetMaskBit("TotalOutlier", false);
					//	qa->SetMaskBit("TerminalOutlier", false);
					//	qa->SetMaskBit("MarginalOutlier", false);
					//	qa->SetMaskBit("SectorLoss", false);
					//	qa->SetMaskBit("LowLiveTime", false);
					//	qa->SetMaskBit("Misc", false);
					//	Keep_event = qa->Pass(run, event_nb);
					// }
				}

				int bad_runs[] = {5610, 5615, 6631, 6757};
				bool Additional_bad_runs = false;
				if (IsData)
				{
					Additional_bad_runs = (std::find(std::begin(bad_runs), std::end(bad_runs), run) != std::end(bad_runs));
				}

				if ((!Keep_event || Additional_bad_runs) && IsData)
					continue;

				////////////////////////////////////////
				/////////QA DB query for charge/////////
				////////////////////////////////////////
				if (IsData)
				{
					qa->AccumulateCharge();
				}

				////////////////////////////////////////////
				// RCDB informations
				////////////////////////////////////////////
				if (IsData && !RCDB_read && !no_RCDB)
				{
					rcdb.readRun(run);
					cout << "Beam energy: " << rcdb.current().beam_energy << endl;
					cout << "Beam current for run " << run << ": " << rcdb.current().beam_current << endl;
					cout << "Beam current requested for run " << run << ": " << rcdb.current().beam_current_request << endl;
					beam_current = rcdb.current().beam_current;
					beam_current_requested = rcdb.current().beam_current_request;
					run_data = run;
					RCDB_read = true;
				}

				///////////////////////////////////////////
				// Get Particles and cut on event topology
				///////////////////////////////////////////
				ev.Set_Particles(PART, IsEE_BG);

				if (ev.recem == 1 && ev.recp == 1)
					Plots.Fill_1D("efficiency", 0, 1);
				if (ev.recep == 1 && ev.recp == 1)
					Plots.Fill_1D("efficiency", 1, 1);
				if (ev.recem == 1 && ev.recep == 1)
					Plots.Fill_1D("efficiency", 2, 1);
				if (ev.recem == 1 && ev.recp == 1 && ev.Proton.status > 2000 && ev.Proton.status < 4000)
					Plots.Fill_1D("efficiency", 3, 1);
				if (ev.recem == 1 && ev.recp == 1 && ev.Proton.status > 4000)
					Plots.Fill_1D("efficiency", 4, 1);
				if (ev.recep == 1 && ev.recp == 1 && ev.Proton.status > 2000 && ev.Proton.status < 4000)
					Plots.Fill_1D("efficiency", 5, 1);
				if (ev.recep == 1 && ev.recp == 1 && ev.Proton.status > 4000)
					Plots.Fill_1D("efficiency", 6, 1);

				if (ev.recem == 1)
					Plots.Fill_1D("efficiency", 7, 1);
				if (ev.recep == 1)
					Plots.Fill_1D("efficiency", 8, 1);
				if (ev.recp == 1)
					Plots.Fill_1D("efficiency", 9, 1);
				if (ev.recp == 1 && ev.Proton.status > 2000 && ev.Proton.status < 4000)
					Plots.Fill_1D("efficiency", 10, 1);
				if (ev.recp == 1 && ev.Proton.status > 4000)
					Plots.Fill_1D("efficiency", 11, 1);

				// Store event topology
				if (ev.recem > 0 && ev.recep > 0 && ev.recp > 0)
				{
					Plots.Fill_1D("multiplicity_e", ev.recem, 1);
					Plots.Fill_1D("multiplicity_ep", ev.recep, 1);
					Plots.Fill_1D("multiplicity_p", ev.recp, 1);

					Plots.Fill_2D("multiplicity_2D_e_ep", ev.recem, ev.recep, 1);
					Plots.Fill_2D("multiplicity_2D_e_p", ev.recem, ev.recp, 1);
					Plots.Fill_2D("multiplicity_2D_ep_p", ev.recep, ev.recp, 1);
				}

				/*if ((ev.topology_FS() != 11 && ev.topology_FS() != 101 && ev.topology_FS() != 110 && ev.topology_FS() != 111) && inclusive_topology)
				{
					continue;
				}*/

				if (!ev.pass_topology_cut() && !inclusive_topology)
				{
					continue;
				}

				///////////////////////////////////////////

				// Number of events after topology cuts
				Plots.Fill_1D("evt_count", 1, 1);

				if (ev.Proton.status > 2000 && ev.Proton.status < 4000)
					Plots.Fill_1D("evt_count", 2, 1);
				if (ev.Proton.status > 4000)
					Plots.Fill_1D("evt_count", 3, 1);

				///////////////////////////////////////////
				// Associate detector responses and do EC cuts
				///////////////////////////////////////////
				ev.Apply_EC_Cuts(CALO);
				ev.Associate_detector_resp(CHE, SCIN);
				ev.Associate_DC_traj(TRAJ);
				ev.Set_Nphe_HTCC();
				///////////////////////////////////////////

				///////////////////////////////////////////
				// TMVA
				///////////////////////////////////////////
				// PositronPID.Evaluate(ev.Positron);
				// ev.Set_Posi_score(PositronPID.score);
				PositronPID.Evaluate_BDT_new(ev.Positron);
				ev.Set_Posi_score(PositronPID.score);

				/*if (!PositronPID.Accept(ev.Positron)) //!!!!!!!!!!!!!!!! Positron cut removed to have consistent lepton ID
					continue;*/

				ElectronPID.Evaluate_BDT_new(ev.Electron);
				ev.Set_Elec_score(ElectronPID.score);
				///////////////////////////////////////////

				// Number of events after positron cuts
				Plots.Fill_1D("evt_count", 4, 1);

				///////////////////////////////////////////
				// Radiative correction
				///////////////////////////////////////////
				ev.Apply_Radiative_Correction(InputParameters.RadCorr);

				///////////////////////////////////////////

				///////////////////////////////////////////
				// Momentum MC correction
				///////////////////////////////////////////
				// cout<<"Before correction "<<ev.Proton.Vector.P()<<endl;  /////////NEED TO BE VALIDATED/////////
				ev.Apply_Energy_loss(EnergyLoss); // Correct proton energy loss
				if (IsData)
				{
					ev.Apply_Momentum_Correction(Momentum_Correction); // Correct momentum of leptons
				}
				// cout<<"After correction "<<ev.Proton.Vector.P()<<endl;
				///////////////////////////////////////////
			}

			nbevent_after_posi++;

			if (!IsHipo)
			{
				input_tree->GetEntry(nbEvent);
				ev.Set_Vectors(*input_Electron, *input_Positron, *input_Proton);
			}

			if (IsData && IsHipo)
			{
				polarization = EVENT.getInt("helicity", 0);
				// polarization = -1. * polarization; ///// TO CHECK //////
			}

			if (!IsData)
			{
				polarization = polarizationGene->Integer(2);
				if (polarization == 0)
					polarization = -1;
			}

			ev.Set_Polarization(polarization);
			ev.Set_Run_Number(run);

			if ((IsHipo) || (!IsHipo))
			{
				nbevent_after_EC++;

				if (HTCCSectorCut)
				{
					if (ev.Positron.SectorCalo(ECAL, PCAL) != ev.Positron.SectorChe(HTCC))
					{
						continue;
					}

					if (ev.Electron.SectorCalo(ECAL, PCAL) != ev.Electron.SectorChe(HTCC))
					{
						continue;
					}
				}

				if (IsHipo)
				{
					///////////////////////////////////////////
					// Momentum Data-driven correction
					///////////////////////////////////////////
					if (false) // if (IsData)
					{
						ev.Apply_Central_Correction(MomCorr, InputParameters.CentralMomCorr, TRACK, TRAJ);
					}

					///////////////////////////////////////////
					// Compute kinematics
					///////////////////////////////////////////
					ev.Get_Kinematics();
					ev.Compute_SF();
				}

				if (!IsHipo)
				{
					ev.Set_Kinematics(input_t, input_Epho, input_qp2, input_M,
									  input_Pt_Frac, input_MMassBeam, input_theta, input_phi, input_s,
									  input_xi, input_L0, input_L);

					ev.Set_SF(input_electron_SF, input_positron_SF);
					ev.Set_Weight(input_weight);
				}

				ev.Get_Polarization_Transfer();

				Plots.Fill_1D("T rec", -ev.t, w);

				// cout<<ev.Positron.TimeChe(HTCC)<<" "<<ev.Electron.TimeChe(HTCC)<<endl;

				outVars["evt_num"] = nbEvent;
				outVars["topology_FS"] = (float)ev.topology_FS();
				outVars["p_p"] = ev.Positron.Vector.P();
				outVars["e_p"] = ev.Electron.Vector.P();
				outVars["prot_p"] = ev.Proton.Vector.P();
				outVars["t"] = ev.t;
				outVars["t_min"] = ev.t_min;
				outVars["MMassBeam"] = ev.MMassBeam;
				outVars["Epho"] = ev.Epho;
				outVars["qp2"] = ev.qp2;
				outVars["M"] = ev.M;
				outVars["xi"] = ev.xi;
				outVars["Pt_Frac"] = ev.Pt_Frac;
				outVars["Q2"] = ev.Q2;
				outVars["theta"] = ev.theta;
				outVars["phi"] = ev.phi;
				outVars["helicity"] = ev.polarization;
				outVars["polaT"] = ev.polaT;
				outVars["positron_SF"] = ev.positron_SF;
				outVars["electron_SF"] = ev.electron_SF;
				outVars["positron_score"] = ev.positron_score;
				outVars["electron_score"] = ev.electron_score;
				outVars["weight"] = ev.w;
				outVars["run"] = ev.run;
				outVars["analysis_stage"] = 0.0;
				outVars["topology"] = (float)ev.topology();
				outVars["positron_Nphe"] = ev.positron_Nphe;
				outVars["electron_Nphe"] = ev.electron_Nphe;
				outVars["positron_HTCCt"] = ev.Positron.TimeChe(HTCC);
				outVars["electron_HTCCt"] = ev.Electron.TimeChe(HTCC);
				outVars["pass_EC_cut"] = ev.pass_EC_cut();
				outVars["positron_HTCC_ECAL_match"] = (ev.Positron.SectorCalo(ECAL, PCAL) == ev.Positron.SectorChe(HTCC)) ? 1. : 0.0;
				outVars["electron_HTCC_ECAL_match"] = (ev.Electron.SectorCalo(ECAL, PCAL) == ev.Electron.SectorChe(HTCC)) ? 1. : 0.0;
				outVars["lead_lep_p"] = (ev.Positron.Vector.P() > ev.Electron.Vector.P()) ? ev.Positron.Vector.P() : ev.Electron.Vector.P();
				outVars["sub_lead_lep_p"] = (ev.Positron.Vector.P() > ev.Electron.Vector.P()) ? ev.Electron.Vector.P() : ev.Positron.Vector.P();
				outVars["lead_lep_theta"] = (ev.Positron.Vector.P() > ev.Electron.Vector.P()) ? ev.Positron.Vector.Theta() : ev.Electron.Vector.Theta();
				outVars["sub_lead_lep_theta"] = (ev.Positron.Vector.P() > ev.Electron.Vector.P()) ? ev.Electron.Vector.Theta() : ev.Positron.Vector.Theta();
				outVars["status_elec"] = ev.Electron.status;
				outVars["status_posi"] = ev.Positron.status;
				outVars["status_prot"] = ev.Proton.status;
				outVars["vx_elec"] = ev.Electron.vertex.x;
				outVars["vy_elec"] = ev.Electron.vertex.y;
				outVars["vz_elec"] = ev.Electron.vertex.z;
				outVars["vx_posi"] = ev.Positron.vertex.x;
				outVars["vy_posi"] = ev.Positron.vertex.y;
				outVars["vz_posi"] = ev.Positron.vertex.z;
				outVars["vx_prot"] = ev.Proton.vertex.z;
				outVars["vy_prot"] = ev.Proton.vertex.z;
				outVars["vz_prot"] = ev.Proton.vertex.z;
				outVars["PCAL_x_elec"] = ev.Electron.X_CALO(PCAL);
				outVars["PCAL_y_elec"] = ev.Electron.Y_CALO(PCAL);
				outVars["PCAL_x_posi"] = ev.Positron.X_CALO(PCAL);
				outVars["PCAL_y_posi"] = ev.Positron.Y_CALO(PCAL);
				outVars["chi2_proton"] = ev.Proton.chi2;
				outVars["Triangular_Cut_elec"] = ((ev.Electron.Energy(ECAL, PCAL) / ev.Electron.Vector.P()) + (ev.Electron.Energy(ECAL, ECIN) / ev.Electron.Vector.P()));
				outVars["Triangular_Cut_posi"] = ((ev.Positron.Energy(ECAL, PCAL) / ev.Positron.Vector.P()) + (ev.Positron.Energy(ECAL, ECIN) / ev.Positron.Vector.P()));
				outVars["CM_gamma_energy"] = CM_gamma_energy(ev.Electron, ev.Positron, ev.Proton, ev.Q2);
				outVars["CM_gamma_energy_2"] = CM_gamma_energy_2(ev.Electron, ev.Positron, ev.Proton);
				outVars["E_k"] = CM_E_k(ev.Electron, ev.Positron, ev.Proton);
				outVars["E_k_2"] = CM_E_k_2(ev.Electron, ev.Positron, ev.Proton, ev.Q2);
				outVars["Q2_true"] = ev.Q2_true;
				outVars["PCAL_sector_elec"] = ev.Electron.SECTOR_CALO(PCAL);
				outVars["PCAL_sector_posi"] = ev.Positron.SECTOR_CALO(PCAL);

				if (!IsData)
				{
					outVars["real_flux_Gen"] = MC_ev.real_flux_Gen;
					outVars["virtual_flux_Gen"] = MC_ev.virtual_flux_Gen;
					outVars["virtual_flux_Frixione_Gen"] = MC_ev.virtual_flux_Frixione_Gen;
					outVars["Epho_Gen"] = MC_ev.Epho_Gen;
					outVars["Q2_Gen"] = MC_ev.Q2_Gen;
					outVars["theta_GJ_Gen"] = MC_ev.theta_GJ_Gen;
					outVars["phi_GJ_Gen"] = MC_ev.phi_GJ_Gen;
				}

				if (PCAL_study)
				{
					outVars["PCAL_U_elec"] = ev.Electron.U_CALO(PCAL);
					outVars["PCAL_V_elec"] = ev.Electron.V_CALO(PCAL);
					outVars["PCAL_W_elec"] = ev.Electron.W_CALO(PCAL);

					outVars["PCAL_U_posi"] = ev.Positron.U_CALO(PCAL);
					outVars["PCAL_V_posi"] = ev.Positron.V_CALO(PCAL);
					outVars["PCAL_W_posi"] = ev.Positron.W_CALO(PCAL);
					// outVars["PCAL_sector_elec"] = ev.Electron.SECTOR_CALO(PCAL);
				}

				if (CALO_study)
				{

					ev.Electron.Get_local_cluster_CALO(PCAL);
					ev.Electron.Get_local_cluster_CALO(ECIN);
					ev.Positron.Get_local_cluster_CALO(PCAL);
					ev.Positron.Get_local_cluster_CALO(ECIN);

					// outVars["PCAL_sector_elec"] = ev.Electron.SECTOR_CALO(PCAL);

					outVars["PCAL_x_elec_rot"] = ev.Electron.cluster_local_PCAL.x;
					outVars["PCAL_y_elec_rot"] = ev.Electron.cluster_local_PCAL.y;
					outVars["PCAL_z_elec_rot"] = ev.Electron.cluster_local_PCAL.z;

					outVars["PCAL_x_posi_rot"] = ev.Positron.cluster_local_PCAL.x;
					outVars["PCAL_y_posi_rot"] = ev.Positron.cluster_local_PCAL.y;
					outVars["PCAL_z_posi_rot"] = ev.Positron.cluster_local_PCAL.z;

					outVars["PCAL_hx_elec_rot"] = ev.Electron.h_cluster_local_PCAL.x;
					outVars["PCAL_hy_elec_rot"] = ev.Electron.h_cluster_local_PCAL.y;
					outVars["PCAL_hz_elec_rot"] = ev.Electron.h_cluster_local_PCAL.z;

					outVars["PCAL_hx_posi_rot"] = ev.Positron.h_cluster_local_PCAL.x;
					outVars["PCAL_hy_posi_rot"] = ev.Positron.h_cluster_local_PCAL.y;
					outVars["PCAL_hz_posi_rot"] = ev.Positron.h_cluster_local_PCAL.z;

					outVars["ECIN_x_elec_rot"] = ev.Electron.cluster_local_ECIN.x;
					outVars["ECIN_y_elec_rot"] = ev.Electron.cluster_local_ECIN.y;
					outVars["ECIN_z_elec_rot"] = ev.Electron.cluster_local_ECIN.z;

					outVars["ECIN_x_posi_rot"] = ev.Positron.cluster_local_ECIN.x;
					outVars["ECIN_y_posi_rot"] = ev.Positron.cluster_local_ECIN.y;
					outVars["ECIN_z_posi_rot"] = ev.Positron.cluster_local_ECIN.z;

					outVars["ECIN_hx_elec_rot"] = ev.Electron.h_cluster_local_ECIN.x;
					outVars["ECIN_hy_elec_rot"] = ev.Electron.h_cluster_local_ECIN.y;
					outVars["ECIN_hz_elec_rot"] = ev.Electron.h_cluster_local_ECIN.z;

					outVars["ECIN_hx_posi_rot"] = ev.Positron.h_cluster_local_ECIN.x;
					outVars["ECIN_hy_posi_rot"] = ev.Positron.h_cluster_local_ECIN.y;
					outVars["ECIN_hz_posi_rot"] = ev.Positron.h_cluster_local_ECIN.z;
				}

				if (Lepton_ID_check)
				{
					outVars["SFPCAL_elec"] = (ev.Electron.Energy(ECAL, PCAL)) / ev.Electron.Vector.P();
					outVars["SFECIN_elec"] = (ev.Electron.Energy(ECAL, ECIN)) / ev.Electron.Vector.P();
					outVars["SFECOUT_elec"] = (ev.Electron.Energy(ECAL, ECOUT)) / ev.Electron.Vector.P();

					outVars["SFPCAL_posi"] = (ev.Positron.Energy(ECAL, PCAL)) / ev.Positron.Vector.P();
					outVars["SFECIN_posi"] = (ev.Positron.Energy(ECAL, ECIN)) / ev.Positron.Vector.P();
					outVars["SFECOUT_posi"] = (ev.Positron.Energy(ECAL, ECOUT)) / ev.Positron.Vector.P();

					outVars["M2PCAL_elec"] = ev.Electron.M2_ECAL(PCAL);
					outVars["M2ECIN_elec"] = ev.Electron.M2_ECAL(ECIN);
					outVars["M2ECOUT_elec"] = ev.Electron.M2_ECAL(ECOUT);

					outVars["M2PCAL_posi"] = ev.Positron.M2_ECAL(PCAL);
					outVars["M2ECIN_posi"] = ev.Positron.M2_ECAL(ECIN);
					outVars["M2ECOUT_posi"] = ev.Positron.M2_ECAL(ECOUT);
				}

				if (DC_Traj_check)
				{
					outVars["DC_R1_elec_x"] = ev.Electron.Trajs[0].x;
					outVars["DC_R1_elec_y"] = ev.Electron.Trajs[0].y;
					outVars["DC_R1_elec_z"] = ev.Electron.Trajs[0].z;
					outVars["DC_R1_posi_x"] = ev.Positron.Trajs[0].x;
					outVars["DC_R1_posi_y"] = ev.Positron.Trajs[0].y;
					outVars["DC_R1_posi_z"] = ev.Positron.Trajs[0].z;

					if (ev.Proton.status < 3000)
					{
						outVars["DC_R1_prot_x"] = ev.Proton.Trajs[0].x;
						outVars["DC_R1_prot_y"] = ev.Proton.Trajs[0].y;
						outVars["DC_R1_prot_z"] = ev.Proton.Trajs[0].z;
					}

					outVars["DC_R2_elec_x"] = ev.Electron.Trajs[1].x;
					outVars["DC_R2_elec_y"] = ev.Electron.Trajs[1].y;
					outVars["DC_R2_elec_z"] = ev.Electron.Trajs[1].z;
					outVars["DC_R2_posi_x"] = ev.Positron.Trajs[1].x;
					outVars["DC_R2_posi_y"] = ev.Positron.Trajs[1].y;
					outVars["DC_R2_posi_z"] = ev.Positron.Trajs[1].z;

					if (ev.Proton.status < 3000)
					{
						outVars["DC_R2_prot_x"] = ev.Proton.Trajs[1].x;
						outVars["DC_R2_prot_y"] = ev.Proton.Trajs[1].y;
						outVars["DC_R2_prot_z"] = ev.Proton.Trajs[1].z;
					}

					outVars["DC_R3_elec_x"] = ev.Electron.Trajs[2].x;
					outVars["DC_R3_elec_y"] = ev.Electron.Trajs[2].y;
					outVars["DC_R3_elec_z"] = ev.Electron.Trajs[2].z;
					outVars["DC_R3_posi_x"] = ev.Positron.Trajs[2].x;
					outVars["DC_R3_posi_y"] = ev.Positron.Trajs[2].y;
					outVars["DC_R3_posi_z"] = ev.Positron.Trajs[2].z;

					if (ev.Proton.status < 3000)
					{
						outVars["DC_R3_prot_x"] = ev.Proton.Trajs[2].x;
						outVars["DC_R3_prot_y"] = ev.Proton.Trajs[2].y;
						outVars["DC_R3_prot_z"] = ev.Proton.Trajs[2].z;
					}
				}

				tree_Electron = ev.Electron.Vector;
				tree_Positron = ev.Positron.Vector;
				tree_Proton = ev.Proton.Vector;
				tree_Missing = ev.vMissing;

				double cutchi2 = 3.;
				// bool chi2cutProton = (vProton.status > 4000 && abs(vProton.chi2 - meanCD) < (cutchi2 * sigmaCD)) || (vProton.status < 4000 && abs(vProton.chi2 - meanFD) < (cutchi2 * sigmaFD));

				if (

					ev.Electron.Vector.P() > InputParameters.ElecMomCut

					&& ev.Positron.Vector.P() > InputParameters.PosiMomCut

					&& abs(ev.Pt_Frac) < InputParameters.PtCut

					&& abs(ev.MMassBeam) < InputParameters.MMassBeamCut

					&& (ev.electron_SF) > InputParameters.ElecMinSF

					&& (ev.positron_SF) > InputParameters.PosiMinSF)
				{
					outVars["analysis_stage"] = 1.0;
					// Number of events after exclusivity cuts
					Plots.Fill_1D("evt_count", 5, 1);

					Plots.Fill_1D("EM1", ev.M, w);

					if (

						ev.Epho < ebeam

						&& ev.Epho > 4.

						&& -ev.t < 0.8

						&& -ev.t > 0.15

						&& ev.qp2 > 1.5 * 1.5

						&& ev.qp2 < 9.

					)
					{

						outVars["analysis_stage"] = 2.0;
						nEventTCS += ev.w;
						denom += (ev.w * ev.w);

						if ((ev.Proton.status / 1000) == 4)
							nCD++;
						if ((ev.Proton.status / 1000) == 2)
							nFD++;

						if ((ev.L / ev.L0) < 0.)
							cout << " L/L0 is negative  " << (ev.L / ev.L0) << "\n";

						float AccValue = Acc_TCS.GetValue(ev.t, ev.Epho, ev.qp2, ev.phi, ev.theta);
						float AccError = Acc_TCS.GetError(ev.t, ev.Epho, ev.qp2, ev.phi, ev.theta);

						ev.Set_Acceptance(AccValue, AccError);
						outVars["acc"] = ev.acc;
						outVars["acc_error"] = ev.acc_error;

						BSA_vs_t.AddEvent(-ev.t, ev.phi, ev.polarization, ev.polaT, AccValue, AccError, ev.w);
						BSA_vs_Xi.AddEvent(ev.xi, ev.phi, ev.polarization, ev.polaT, AccValue, AccError, ev.w);
						BSA_vs_M.AddEvent(sqrt(ev.qp2), ev.phi, ev.polarization, ev.polaT, AccValue, AccError, ev.w);

						RRatio_vs_t.AddEvent(-ev.t, ev.phi, ev.theta, ev.L, ev.L0, AccValue, AccError, ev.w);
						RRatio_vs_Xi.AddEvent(ev.xi, ev.phi, ev.theta, ev.L, ev.L0, AccValue, AccError, ev.w);

						int binningt = Acc_TCS.Get_t_bin(ev.t);
						int binningMass = Acc_TCS.Get_M_bin(ev.qp2);
						int binningEg = Acc_TCS.Get_Eg_bin(ev.Epho);

						float CorrVolume1 = Bin_Volume_Corr.volume1[binningt + Acc_TCS.bin_t * binningMass + Acc_TCS.bin_t * Acc_TCS.bin_M * binningEg];
						float CorrVolume2 = Bin_Volume_Corr.volume2[binningt + Acc_TCS.bin_t * binningMass + Acc_TCS.bin_t * Acc_TCS.bin_M * binningEg];
						AFB_vs_t.AddEvent(-ev.t, ev.phi, ev.theta, AccValue, AccError, ev.w, CorrVolume1, CorrVolume2);
						AFB_vs_t1.AddEvent(-ev.t, ev.phi, ev.theta, AccValue, AccError, ev.w, CorrVolume1, CorrVolume2);
						AFB_vs_M.AddEvent(sqrt(ev.qp2), ev.phi, ev.theta, AccValue, AccError, ev.w, CorrVolume1, CorrVolume2);
						AFB_vs_Eg.AddEvent(ev.Epho, ev.phi, ev.theta, AccValue, AccError, ev.w, CorrVolume1, CorrVolume2);
						if (ev.qp2 < 4.)
							AFB_vs_t2.AddEvent(-ev.t, ev.phi, ev.theta, AccValue, AccError, ev.w, CorrVolume1, CorrVolume2);

						if (AccValue <= 0.0)
							cout << "Acc Value bellow zero : problem"
								 << "\n";
					}
				}

				outT->Fill();
			}
		}
	}

	/*AFB_vs_t.Calculate();
	AFB_vs_t1.Calculate();
	AFB_vs_M.Calculate();
	AFB_vs_Eg.Calculate();
	AFB_vs_t2.Calculate();

	BSA_vs_t.Calculate(InputParameters.factorPola);
	BSA_vs_Xi.Calculate(InputParameters.factorPola);
	BSA_vs_M.Calculate(InputParameters.factorPola);

	RRatio_vs_t.Calculate();
	RRatio_vs_Xi.Calculate();*/

	///////////////////////////////
	// Draw Plots here
	///////////////////////////////
	Plots.Draw_All_1D();
	Plots.Draw_All_2D();

	outFile->cd();
	outT->Write();
	outFile->Write();
	outFile->Close();

	cout << "Tree written" << endl;

	///////////////////////////////
	// Output charge from QA in txt
	///////////////////////////////
	if (IsData)
	{
		std::ofstream output_txt_charge(Form("outputTCS_charge_%s.txt", output_file.Data()));
		if (!output_txt_charge.is_open())
		{
			std::cerr << "Unable to open the charge file." << std::endl;
		}
		output_txt_charge << "Charge" << std::endl;
		output_txt_charge << (qa->GetAccumulatedCharge()) << std::endl;
		output_txt_charge << "Current" << std::endl;
		output_txt_charge << beam_current << std::endl;
		output_txt_charge << "Current requested" << std::endl;
		output_txt_charge << beam_current_requested << std::endl;
		output_txt_charge << "run" << std::endl;
		output_txt_charge << run_data << std::endl;
	}

	cout << "nb of file " << nbf << "\n";
	cout << "nb of events " << nbEvent << "\n";
	cout << "real number od equivalent tcs " << (nEventTCS * nEventTCS) / denom << "\n";
	cout << "nb events after positron cut " << nbevent_after_posi << "\n";
	cout << "nb events after EC cut " << nbevent_after_EC << "\n";
	cout << "nb event TCS " << nEventTCS << "\n";
	cout << "nb event CD " << nCD << "\n";
	cout << "nb event FD " << nFD << "\n";
	cout << "AfterCuts " << AfterCuts << "\n";
	cout << "nb JPSI " << nbJPSI << "\n";

	// gROOT->ProcessLine(".q");

	time(&end); // note time after execution

	double difference = difftime(end, begin);
	printf("All this work done in only %.2lf seconds. Congratulations !\n", difference);

	gApplication->Terminate();

	return 0;
}
