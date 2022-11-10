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
#include "bib/TCSRunSelector.h"

#include "reader.h"

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

	bool IsData = true;
	bool IsHipo = true;

	bool IsTCSGen = false;
	bool IsGrape = true;
	bool IsJPsi = false;
	bool Weighted_simu = false;

	bool HTCCSectorCut = false;

	bool RGA_Fall2018 = false;

	Int_t argc = gApplication->Argc();
	char **argv = gApplication->Argv();
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

	//////////////////////////////////////////////
	// Acceptance setup
	//////////////////////////////////////////////

	Acceptance Acc_TCS(TString(argv[argc - 3]), 4, 3, 3, 36, 13);
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
	TString output_file = (TString)(argv[argc - 1]);
	TFile *outFile = new TFile(Form("outputTCS_%s.root", output_file.Data()), "recreate");
	// TFile *outFile = new TFile("outputTCS_"+output_file+".root", "recreate");
	TTree *outT = new TTree("tree", "tree");

	TLorentzVector tree_Electron, tree_Positron, tree_Proton, tree_Missing;
	outT->Branch("Electron", "TLorentzVector", &tree_Electron);
	outT->Branch("Positron", "TLorentzVector", &tree_Positron);
	outT->Branch("Proton", "TLorentzVector", &tree_Proton);
	outT->Branch("Missing", "TLorentzVector", &tree_Missing);

	int trigger_bit;
	outT->Branch("trigger_bit", &trigger_bit, "trigger_bit/I");

	TString fvars[] = {
		"t", "MMassBeam", "Epho", "qp2", "M", "xi", "s", "L", "L0", "Pt_Frac", "Q2", "theta", "phi", "positron_SF", "electron_SF", "positron_score", "electron_score",
		"weight", "acc", "acc_error", "real_flux", "virtual_flux", "run", "analysis_stage", "topology",
		"positron_Nphe", "electron_Nphe", "positron_HTCCt", "electron_HTCCt", "positron_HTCC_ECAL_match", "electron_HTCC_ECAL_match"};

	std::map<TString, Float_t>
		outVars;
	for (size_t i = 0; i < sizeof(fvars) / sizeof(TString); i++)
	{
		outVars[fvars[i]] = 0.;
		ADDVAR(&(outVars[fvars[i]]), fvars[i], "/F", outT);
	}

	TString fvars_Gen[] = {
		"t_Gen", "MMassBeam_Gen", "Epho_Gen", "qp2_Gen", "M_Gen_1", "M_Gen_2", "Pt_Frac_Gen", "Q2_Gen", "theta_Gen", "phi_Gen"};

	if (IsGrape || IsTCSGen)
	{
		for (size_t i = 0; i < sizeof(fvars_Gen) / sizeof(TString); i++)
		{
			cout<<"here\n";
			outVars[fvars_Gen[i]] = 0.;
			ADDVAR(&(outVars[fvars_Gen[i]]), fvars_Gen[i], "/F", outT);
		}
	}
	///////////////////////////////////////////

	int nbtc = 0;
	int nbJPSI = 0;

	///////////////////////////////////////////
	// TMVA PID for Positron
	///////////////////////////////////////////
	PositronIdentification PositronPID("MLP method", "TMVAClassification_MLP6D.weights.xml", InputParameters.MLPscoreCut, InputParameters.MLPMomCut);
	PositronPID.InitializePositronIdentification();

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

	int seed = atoi(argv[5]);
	TRandom *polarizationGene = new TRandom(seed);
	TRandom3 *ChoiceEvent = new TRandom3(seed);

	int corrrad = 0;
	int nbEvent = 0;

	////////////////////////////////////////////
	// Get file name
	////////////////////////////////////////////
	for (Int_t i = 3; i < (argc - 3); i++)
	{
		if (TString(argv[i]).Contains("MC"))
		{
			IsData = false;
		}

		if (TString(argv[i]).Contains(".hipo"))
		{
			nbf++;
			nameFiles = TString(argv[i]);
		}

		if (TString(argv[3]).Contains(".root"))
		{
			IsHipo = false;
			nameFiles = TString(argv[i]);
		}

		if (IsTCSGen || IsGrape || IsJPsi)
		{
			IsData = false;
		}

		////////////////////////////////////////////
		cout << "////////////////////////////////////////////" << "\n";
		if (IsData)
			cout << "Running on Data" << "\n";
		else if (IsTCSGen)
			cout << "Running on TCSGen Simulation" << "\n";
		else if (IsGrape)
			cout << "Running on Grape Simulation" << "\n";
		else if (IsJPsi)
			cout << "Running on JPsi Simulation" << "\n";
		cout << "////////////////////////////////////////////" << "\n";
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
			factory.show();
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

		while (((reader.next() && IsHipo) || (nbEvent < nentries && !IsHipo)) /*&& nbEvent < 10000*/)
		{

			nbEvent++;
			if (nbEvent % 30000 == 0)
			{
				time(&intermediate);
				double intermediate_time = difftime(intermediate, begin);

				cout << nbEvent << " events processed in " << intermediate_time << "s" << "\n";
			}
			// Number of total event
			Plots.Fill_1D("evt_count", 0, 1);
			// cout << "event" << endl;
			double w = 1; // MCfluxBH[0]*MCpsfBH[0]*MCcsBH[0];
			int polarization;
			Event ev;
			MCEvent MC_ev;

			int run = 0;
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

				if (PART.getSize() < 1)
					continue;

				run = RUN.getInt("run", 0);
				trigger_bit = RUN.getLong("trigger", 0);
				int np_input = PART.getRows();
				ev.Set_nb_part(np_input);
				ev.Set_trigger_bit(trigger_bit);

				MC_ev.Set_MC_Particles(MCEVENT, MCPART, IsGrape);
				MC_ev.Get_Kinematics();

				if (!IsData)
					w = MC_ev.w;

				ev.Set_Weight(w);

				///////////////////////////////////////////
				// Filter good runs
				///////////////////////////////////////////
				if (!Run_Selector.Is_Good_Run(run) && IsData && RGA_Fall2018)
					continue;

				///////////////////////////////////////////
				// Get Particles and cut on event topology
				///////////////////////////////////////////
				ev.Set_Particles(PART);
				if (!ev.pass_topology_cut())
				{
					continue;
				}
				// cout<<"event"<<endl;
				// cout<<ev.topology()<<endl;
				// PART.show();
				///////////////////////////////////////////

				// Number of events after topology cuts
				Plots.Fill_1D("evt_count", 1, 1);

				///////////////////////////////////////////
				// Associate detector responses and do EC cuts
				///////////////////////////////////////////
				ev.Apply_EC_Cuts(CALO);
				ev.Associate_detector_resp(CHE, SCIN);
				ev.Set_Nphe_HTCC();
				///////////////////////////////////////////

				///////////////////////////////////////////
				// TMVA
				///////////////////////////////////////////
				PositronPID.Evaluate(ev.Positron);
				ev.Set_Posi_score(PositronPID.score);
				// cout<<"event"<<endl;
				// cout<<PositronPID.score<<endl;
				//  cout<<PositronPID.Accept(ev.Positron)<<endl;

				if (!PositronPID.Accept(ev.Positron))
					continue;

				PositronPID.Evaluate(ev.Electron);
				ev.Set_Elec_score(PositronPID.score);
				///////////////////////////////////////////

				// Number of events after positron cuts
				Plots.Fill_1D("evt_count", 2, 1);

				///////////////////////////////////////////
				// Radiative correction
				///////////////////////////////////////////
				ev.Apply_Radiative_Correction(InputParameters.RadCorr);
				ev.Compute_SF();
				///////////////////////////////////////////

				///////////////////////////////////////////
				// Momentum MC correction
				///////////////////////////////////////////
				ev.Apply_MC_Correction(MomCorr);
				///////////////////////////////////////////
			}
			// cout << nbevent_after_posi << endl;
			/*if (pass_positron_cut == false)
				cout << pass_positron_cut << endl;
*/

			nbevent_after_posi++;
			// cout << nbevent_after_posi << endl;

			if (!IsHipo)
			{
				// cout<<"here"<<endl;
				input_tree->GetEntry(nbEvent);
				// cout<<"here"<<endl;
				ev.Set_Vectors(*input_Electron, *input_Positron, *input_Proton);
			}

			if (IsData && IsHipo)
			{
				polarization = EVENT.getInt("helicity", 0);
				polarization = -1. * polarization; ///// TO CHECK //////
			}

			if (!IsData)
			{
				polarization = polarizationGene->Integer(2);
				if (polarization == 0)
					polarization = -1;
			}

			ev.Set_Polarization(polarization);
			ev.Set_Run_Number(run);

			if ((ev.pass_EC_cut() && IsHipo) || (!IsHipo))
			{
				nbevent_after_EC++;

				if (HTCCSectorCut)
				{
					if (ev.Positron.SectorCalo(ECAL, PCAL) != ev.Positron.SectorChe(HTCC))
					{
						// cout<<"mismatch posi"<<endl;
						continue;
					}

					if (ev.Electron.SectorCalo(ECAL, PCAL) != ev.Electron.SectorChe(HTCC))
					{
						// cout<<"mismatch elec"<<endl;
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

				outVars["p_p"] = ev.Positron.Vector.P();
				outVars["e_p"] = ev.Electron.Vector.P();
				outVars["prot_p"] = ev.Proton.Vector.P();
				outVars["t"] = ev.t;
				outVars["MMassBeam"] = ev.MMassBeam;
				outVars["Epho"] = ev.Epho;
				outVars["qp2"] = ev.qp2;
				outVars["M"] = ev.M;
				outVars["xi"] = ev.xi;
				outVars["Pt_Frac"] = ev.Pt_Frac;
				outVars["Q2"] = ev.Q2;
				outVars["theta"] = ev.theta;
				outVars["phi"] = ev.phi;
				outVars["positron_SF"] = ev.positron_SF;
				outVars["electron_SF"] = ev.electron_SF;
				outVars["positron_score"] = ev.positron_score;
				outVars["electron_score"] = ev.electron_score;
				outVars["weight"] = ev.w;
				outVars["real_flux"] = ev.real_flux;
				outVars["virtual_flux"] = ev.virtual_flux;
				outVars["run"] = ev.run;
				outVars["analysis_stage"] = 0.0;
				outVars["topology"] = (float)ev.topology();
				outVars["positron_Nphe"] = ev.positron_Nphe;
				outVars["electron_Nphe"] = ev.electron_Nphe;
				outVars["positron_HTCCt"] = ev.Positron.TimeChe(HTCC);
				outVars["electron_HTCCt"] = ev.Electron.TimeChe(HTCC);
				outVars["positron_HTCC_ECAL_match"] = (ev.Positron.SectorCalo(ECAL, PCAL) == ev.Positron.SectorChe(HTCC)) ? 1. : 0.0;
				outVars["electron_HTCC_ECAL_match"] = (ev.Electron.SectorCalo(ECAL, PCAL) == ev.Electron.SectorChe(HTCC)) ? 1. : 0.0;

				if (IsGrape || IsTCSGen)
				{
					outVars["t_Gen"] = MC_ev.t_Gen;
					outVars["MMassBeam_Gen"] = MC_ev.MMassBeam_Gen;
					outVars["Epho_Gen"] = MC_ev.Epho_Gen;
					outVars["qp2_Gen"] = MC_ev.qp2_Gen;
					outVars["M_Gen_1"] = MC_ev.M_Gen_1;
					outVars["M_Gen_2"] = MC_ev.M_Gen_2;
					outVars["Pt_Frac_Gen"] = MC_ev.Pt_Frac_Gen;
					outVars["Q2_Gen"] = MC_ev.Q2_Gen;
					outVars["theta_Gen"] = MC_ev.theta_Gen;
					outVars["phi_Gen"] = MC_ev.phi_Gen;
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
					Plots.Fill_1D("evt_count", 3, 1);

					Plots.Fill_1D("EM1", ev.M, w);

					if (

						ev.Epho < 10.6

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
							cout << "Acc Value bellow zero : problem" << "\n";
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

	outT->Write();
	outFile->Write();
	outFile->Close();

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
