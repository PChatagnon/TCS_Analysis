#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "THStack.h"
#include "bib_CS_extraction/Fit_Function_Class.h"
#include "bib_CS_extraction/Table_Class.h"
#include "bib_CS_extraction/Sample_Class.h"
#include "bib_CS_extraction/Analysis_Class.h"
#include "bib_CS_extraction/Utils.h"
#include "../bib/InputParser.h"
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

int CS_Extraction_combine()
{

	time_t begin, filtering_time, fitting_time, MC_time, end; 
	time(&begin);

	gROOT->SetBatch(kTRUE);
	gStyle->SetPaintTextFormat("4.1f");
	gStyle->SetOptStat(1);
	gStyle->SetPalette(55);
	gStyle->SetLabelSize(.04, "xyz");
	gStyle->SetTitleSize(.04, "xyz");
	gStyle->SetTitleSize(.07, "t");
	gStyle->SetFrameLineWidth(0);
	gStyle->SetLineWidth(1);
	gStyle->SetHistLineWidth(1);
	gStyle->SetMarkerStyle(13);
	gStyle->SetTitleW(0.8); // per cent of the pad width
	gStyle->SetTitleH(0.1); // per cent of the pad height

	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
	ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);


	///////////////////////////////////////////////////////
	////////////// Instantiate analysis class /////////////
	///////////////////////////////////////////////////////
	Analysis JPsi_CS_analysis;

	///////////////////////////////////////////////////////
	//////// Read input configuration file ////////////////
	///////////////////////////////////////////////////////
	Int_t argc = gApplication->Argc();
	char **argv = gApplication->Argv();
	Input input(argc, argv);

	string parameters_file;
	if (input.cmdOptionExists("-param"))
	{
		parameters_file = input.getCmdOption("-param");
		JPsi_CS_analysis.Set_Parameters(parameters_file);
	}

	///////////////////////////////////////////////////////
	////////////////   Setup the analysis  ////////////////
	///////////////////////////////////////////////////////
	JPsi_CS_analysis.Set_Sample_to_RGA();
	JPsi_CS_analysis.Set_Binning_RGA();

	// JPsi_CS_analysis.Set_Sample_to_RGB();
	// JPsi_CS_analysis.Set_Binning_RGB();

	JPsi_CS_analysis.Setup_Histo();
	JPsi_CS_analysis.Setup_Latex_Table();
	///////////////////////////////////////////////////////

	////////////////////////////////////
	//////// Filter the samples ////////
	////////////////////////////////////
	JPsi_CS_analysis.Analysis_Sample.Filter_Sample_Trees(JPsi_CS_analysis.kinematic_cut, JPsi_CS_analysis.kinematic_cut_BG, JPsi_CS_analysis.exclusivity_cut);
	JPsi_CS_analysis.Analysis_Sample.Filter_Data_Tree(JPsi_CS_analysis.kinematic_cut, JPsi_CS_analysis.data_cut, JPsi_CS_analysis.exclusivity_cut);
	time(&filtering_time);
	////////////////////////////////////

	/////////////////////////////////////////
	//////// Cross section calculation //////
	/////////////////////////////////////////
	JPsi_CS_analysis.Process_REC();
	time(&fitting_time);
	JPsi_CS_analysis.Process_MC();
	time(&MC_time);
	JPsi_CS_analysis.Finalize_int_CS_Calculation();

	JPsi_CS_analysis.Save_Histo_to_root();
	JPsi_CS_analysis.Latex_Table.Format();
	JPsi_CS_analysis.Latex_Table.Save();
	/////////////////////////////////////////


	///////////////////////////////////
	///// Final message and close /////
	///////////////////////////////////
	time(&end);
	double total_t = difftime(end, begin);
	double filtering_t = difftime(filtering_time, begin);
	double fitting_t = difftime(fitting_time, filtering_time);
	double MC_t = difftime(MC_time, fitting_time);
	double finalizing_t = difftime(end, MC_time);

	cout << " Time summary: " << endl;
	printf("Total time: %.2lf seconds\n", total_t);
	printf("Filtering time: %.2lf seconds\n", filtering_t);
	printf("Fitting time: %.2lf seconds\n", fitting_t);
	printf("MC processing time: %.2lf seconds\n", MC_t);
	printf("Final plotting time: %.2lf seconds\n", finalizing_t);
	///////////////////////////////////

	gApplication->Terminate();

	return 0;
}
