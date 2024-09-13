#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "THStack.h"
#include <iostream>
#include <fstream>
#include "../CS_Extraction/bib_CS_extraction/Table_Class.h"
using namespace std;

void resetGraphYValues(TGraphAsymmErrors *graph)
{
	if (!graph)
		return; // Check if the graph pointer is valid

	// Loop over all points in the graph
	for (int i = 0; i < graph->GetN(); i++)
	{
		double graph_x, graph_y;
		graph->GetPoint(i, graph_x, graph_y); // Get the current point (x, y)
		graph->SetPoint(i, graph_x, 0.0);	  // Set the y value to 0.0
		graph->SetPointError(i, 0.1, 0.1, 0.0, 0.0);
	}
}

void updateGraphErrors(TGraphAsymmErrors *graph_base, TGraphAsymmErrors *graph_down, TGraphAsymmErrors *graph_up)
{
	if (!graph_base || !graph_down || !graph_up)
		return; // Check if the pointers are valid

	int nPoints = graph_base->GetN();
	int nPointsNewErrors = graph_down->GetN();

	// Ensure both graphs have the same number of points
	if (nPoints != nPointsNewErrors)
	{
		std::cerr << "Error: The number of points in both graphs does not match!" << std::endl;
		return;
	}

	// Loop over all points in the graph
	for (int i = 0; i < nPoints; i++)
	{
		double ErrYLow, ErrYHigh, newErrYLow_down, newErrYLow_up, newErrYLow, newErrYHigh;

		ErrYLow = graph_base->GetErrorYlow(i);
		ErrYHigh = graph_base->GetErrorYhigh(i);

		newErrYLow_down = -1.0 * graph_down->GetErrorYlow(i);
		cout << "Down " << newErrYLow_down << endl;

		newErrYLow_up = -1.0 * graph_up->GetErrorYlow(i);
		cout << "Up " << newErrYLow_up << endl;

		newErrYLow = std::min({newErrYLow_down, newErrYLow_up, 0.0});
		newErrYHigh = std::max({newErrYLow_down, newErrYLow_up, 0.0});
		cout << "Here " << newErrYLow << " " << newErrYHigh << endl;

		cout << "Initial errors " << ErrYLow << " " << ErrYHigh << endl;
		cout << "" << endl;

		// Update the errors in the main graph with the new errors
		graph_base->SetPointError(i, 0.1, 0.1, sqrt(ErrYLow * ErrYLow + newErrYLow * newErrYLow), sqrt(ErrYHigh * ErrYHigh + newErrYHigh * newErrYHigh));
	}
}

void printGraphContent(TGraphAsymmErrors* graph) {
    if (!graph) {
        std::cerr << "Error: Invalid TGraphAsymmErrors pointer!" << std::endl;
        return; 
    }

    int nPoints = graph->GetN(); 
    std::cout << "Low Errors: ";
    
    for (int i = 0; i < nPoints; ++i) {
        double eyl = graph->GetErrorYlow(i);
        std::cout << eyl << ", ";
    }

    std::cout << std::endl; // End line for X errors

    std::cout << "High Errors: ";
    
    // Loop over all points to collect and print Y errors
    for (int i = 0; i < nPoints; ++i) {
        double eyh = graph->GetErrorYhigh(i);
		std::cout << eyh << ", ";
    }
    std::cout << std::endl; // End line for Y errors
}

int Systematic_Variation_JPsi()
{
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

	std::vector<vector<TString>> Systematics_array{

		/*{"Name","Title", "name_base", "File_Base", "name_down", "File_Down", "name_up", "File_Up"},*/

		{"Q2_cut", "Q2 cut",
		 "0.5 GeV^{2}", "../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root", ";1",
		 "0.2 GeV^{2}", "../CS_Extraction/Results_CS/CS_Q2_Systematic/CS_Extraction_combine_Q2_02_CS_graph.root", ";1", "1.0",
		 "0.8 GeV^{2}", "../CS_Extraction/Results_CS/CS_Q2_Systematic/CS_Extraction_combine_Q2_08_CS_graph.root", ";1", "1.0",
		 "70", "-70", "E_{#gamma} [GeV]", "log", "int"},

		{"MM_Cut", "Missing mass cut",
		 "0.4 GeV^{2}", "../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root", ";1",
		 "0.2 GeV^{2}", "../CS_Extraction/Results_CS/CS_MM_Systematic/CS_Extraction_combine_MM_02_CS_graph.root", ";1", "1.",
		 "0.8 GeV^{2}", "../CS_Extraction/Results_CS/CS_MM_Systematic/CS_Extraction_combine_MM_08_CS_graph.root", ";1", "1.",
		 "10", "-10", "E_{#gamma} [GeV]", "log", "int"},

		{"Fit", "Fit function",
		 "Exp. BG", "../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root", ";1",
		 "Fit 1", "../CS_Extraction/Results_CS/CS_Fit_Systematic/CS_Extraction_combine_Fit_1_CS_graph.root", ";1", "1.",
		 "Fit 4", "../CS_Extraction/Results_CS/CS_Fit_Systematic/CS_Extraction_combine_Fit_4_CS_graph.root", ";1", "1.",
		 "50", "-50", "E_{#gamma} [GeV]", "log", "int"},

		{"AI_PID", "AI PID",
		 "Default", "../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root", ";1",
		 "AI PID 1", "../CS_Extraction/Results_CS/CS_AI_PID_Systematic/CS_Extraction_combine_AI_PID_Systematic_1_CS_graph.root", ";1", "1.",
		 "AI PID 2", "../CS_Extraction/Results_CS/CS_AI_PID_Systematic/CS_Extraction_combine_AI_PID_Systematic_2_CS_graph.root", ";1", "1.",
		 "50", "-50", "E_{#gamma} [GeV]", "log", "int"},

		{"Lep_Mom", "Lepton Momentum cut",
		 "Lepton mom. 1.7", "../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root", ";1",
		 "Lepton mom. 1.5", "../CS_Extraction/Results_CS/CS_Mom_lepton_Systematic/CS_Extraction_combine_Mom_lepton_1_5_CS_graph.root", ";1", "1.",
		 "Lepton mom. 1.9", "../CS_Extraction/Results_CS/CS_Mom_lepton_Systematic/CS_Extraction_combine_Mom_lepton_1_9_CS_graph.root", ";1", "1.",
		 "50", "-50", "E_{#gamma} [GeV]", "log", "int"},

		{"Proton_PID", "Proton PID",
		 "No cut", "../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root", ";1",
		 "2-sigma", "../CS_Extraction/Results_CS/CS_Proton_PID_Systematic/CS_Extraction_combine_nominal_Proton_PID_1_CS_graph.root", ";1", "1.",
		 "3-sigma", "../CS_Extraction/Results_CS/CS_Proton_PID_Systematic/CS_Extraction_combine_nominal_Proton_PID_2_CS_graph.root", ";1", "1.",
		 "50", "-50", "E_{#gamma} [GeV]", "log", "int"},

		/*

		{"Q2_cut", "Q2 cut",
		 "0.5 GeV^{2}", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_with_int_SLER_CS_graph.root", "Graph",
		 "0.2 GeV^{2}", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_Q2_02_CS_graph.root", "Graph","1.25",
		 "0.8 GeV^{2}", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_Q2_10_CS_graph.root", "Graph","0.92",
		 "50", "-50", "E_{#gamma} [GeV]", "no_log"},

		{"MM_Cut", "Missing mass cut",
		 "0.4 GeV^{2}", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_with_int_SLER_CS_graph.root", "Graph",
		 "0.2 GeV^{2}", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_MM_08_CS_graph.root", "Graph","1.",
		 "0.8 GeV^{2}", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_MM_08_CS_graph.root", "Graph","1.",
		 "10", "-10", "E_{#gamma} [GeV]", "no_log"},

		{"Rad_Corr", "Rad. Corr.",
		 "With BH rad. corr.", "../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root", ";1",
		 "With VM rad. corr", "../CS_Extraction/CS_Extraction_Gauss_exp_1_CS_graph.root", ";1","1.",
		 "With VM rad. corr", "../CS_Extraction/CS_Extraction_Gauss_exp_1_CS_graph.root", ";1","1.",
		 "5", "-15", "E_{#gamma} [GeV]", "no_log"},

		{"Fit", "Fit function",
		 "Exp. BG", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_with_int_SLER_CS_graph.root", "Graph",
		 "Pol. 2 BG", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_with_int_SLER_pol_BG_CS_graph.root", "Graph","1.",
		 "Pol. 2 BG", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_with_int_SLER_pol_BG_CS_graph.root", "Graph","1.",
		 "50", "-50", "E_{#gamma} [GeV]", "no_log"},

		{"Flux", "Flux fonction",
		 "Q_{min}<0.02", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_with_int_SLER_CS_graph.root", "Graph",
		 "Frixione (no cut)", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_Frixione_Flux_1_CS_graph.root", "Graph","1.",
		 "Frixione (no cut)", "../CS_Extraction/CS_Extraction_combine_good_error_Gaussian_Frixione_Flux_1_CS_graph.root", "Graph","1.",
		 "50", "-50", "E_{#gamma} [GeV]", "no_log"},

		 */

		/*{"Rad_Corr", "Rad. Corr. - t-dependence ",
		 "With rad. corr.", "../CS_Extraction/t_cross_section/CS_Extraction_combine_t_05_rad_new_9.28_10.36_CS_graph.root", "Graph",
		 "No Rad.", "../CS_Extraction/t_cross_section/CS_Extraction_combine_t_05_rad_new_no_rad9.28_10.36_CS_graph.root", "Graph",
		 "No Rad.", "../CS_Extraction/t_cross_section/CS_Extraction_combine_t_05_rad_new_no_rad9.28_10.36_CS_graph.root", "Graph",
		 "10", "-10", "-t [GeV^{2}]", "log"},*/

		{"M2_cut", "MM2 cut - bin 1 - t-dependence",
		 "0.5", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_8.20_9.28_CS_graph.root", ";1",
		 "0.2", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/MM_02_config_bin1_8.20_9.28_CS_graph.root", ";1", "1.",
		 "0.8", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/MM_08_config_bin1_8.20_9.28_CS_graph.root", ";1", "1.",
		 "50", "-50", "-t [GeV^{2}]", "log", "bin1"},

		{"M2_cut", "MM2 cut - bin 2 - t-dependence",
		 "0.5", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_9.28_10.00_CS_graph.root", ";1",
		 "0.2", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/MM_02_config_bin2_9.28_10.00_CS_graph.root", ";1", "1.",
		 "0.8", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/MM_08_config_bin2_9.28_10.00_CS_graph.root", ";1", "1.",
		 "50", "-50", "-t [GeV^{2}]", "log", "bin2"},

		{"M2_cut", "MM2 cut - bin 3 - t-dependence",
		 "0.5", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_10.00_10.60_CS_graph.root", ";1",
		 "0.2", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/MM_02_config_bin3_10.00_10.60_CS_graph.root", ";1", "1.",
		 "0.8", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/MM_08_config_bin3_10.00_10.60_CS_graph.root", ";1", "1.",
		 "50", "-50", "-t [GeV^{2}]", "log", "bin3"},

		{"Q2_cut", "Q2 cut - bin 1 - t-dependence",
		 "0.5", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_8.20_9.28_CS_graph.root", ";1",
		 "0.2", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/Q2_02_config_bin1_8.20_9.28_CS_graph.root", ";1", "1.25",
		 "0.8", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/Q2_08_config_bin1_8.20_9.28_CS_graph.root", ";1", "0.92",
		 "50", "-50", "-t [GeV^{2}]", "log", "bin1"},

		{"Q2_cut", "Q2 cut - bin 2 - t-dependence",
		 "0.5", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_9.28_10.00_CS_graph.root", ";1",
		 "0.2", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/Q2_02_config_bin2_9.28_10.00_CS_graph.root", ";1", "1.25",
		 "0.8", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/Q2_08_config_bin2_9.28_10.00_CS_graph.root", ";1", "0.92",
		 "50", "-50", "-t [GeV^{2}]", "log", "bin2"},

		{"Q2_cut", "Q2 cut - bin 3 - t-dependence",
		 "0.5", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_10.00_10.60_CS_graph.root", ";1",
		 "0.2", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/Q2_02_config_bin3_10.00_10.60_CS_graph.root", ";1", "1.25",
		 "0.8", "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/Q2_08_config_bin3_10.00_10.60_CS_graph.root", ";1", "0.92",
		 "50", "-50", "-t [GeV^{2}]", "log", "bin3"},

	};

	// Store the combinaison of systematics here
	TFile *total_graph_int_file = new TFile("../CS_Extraction/Results_CS/CS_Nominal/CS_Extraction_combine_nominal_CS_graph.root");
	auto total_graph_int = (TGraphAsymmErrors *)total_graph_int_file->Get(";1");
	resetGraphYValues(total_graph_int);

	TFile *total_graph_bin1_file = new TFile("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_8.20_9.28_CS_graph.root");
	auto total_graph_bin1 = (TGraphAsymmErrors *)total_graph_bin1_file->Get(";1");
	resetGraphYValues(total_graph_bin1);

	TFile *total_graph_bin2_file = new TFile("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_9.28_10.00_CS_graph.root");
	auto total_graph_bin2 = (TGraphAsymmErrors *)total_graph_bin2_file->Get(";1");
	resetGraphYValues(total_graph_bin2);

	TFile *total_graph_bin3_file = new TFile("/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction/t_cross_section/Results_diff_CS/CS_Nominal/CS_Extraction_combine_t_05_rad_new_2_10.00_10.60_CS_graph.root");
	auto total_graph_bin3 = (TGraphAsymmErrors *)total_graph_bin3_file->Get(";1");
	resetGraphYValues(total_graph_bin3);

	for (int i = 0; i < Systematics_array.size(); i++)
	{

		//////////////////////////////////////////////////
		// Set options for each label
		TString name_sys = Systematics_array[i][0];
		TString title_sys = Systematics_array[i][1];
		TString name_base = Systematics_array[i][2];
		TString file_base = Systematics_array[i][3];
		TString name_file_base = Systematics_array[i][4];
		TString name_down = Systematics_array[i][5];
		TString file_down = Systematics_array[i][6];
		TString name_file_down = Systematics_array[i][7];
		double weight_file_down = std::stod(Systematics_array[i][8].Data());
		TString name_up = Systematics_array[i][9];
		TString file_up = Systematics_array[i][10];
		TString name_file_up = Systematics_array[i][11];
		double weight_file_up = std::stod(Systematics_array[i][12].Data());

		double max_y_syst = std::stod(Systematics_array[i][13].Data());
		double min_y_syst = std::stod(Systematics_array[i][14].Data());
		TString var_label = Systematics_array[i][15];
		TString log_label = Systematics_array[i][16];

		TString combinaison_label = Systematics_array[i][17];

		//////////////////////////////////////////////////

		cout << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << "Doing " << title_sys << " systematic plot" << endl;
		cout << "//////////////////////////////////////////////////" << endl;
		cout << endl;

		////////////////////////Get Graphs //////////////////////////////
		TFile *base_file = new TFile(file_base);
		auto base_graph = (TGraphAsymmErrors *)base_file->Get(name_file_base);

		TFile *down_file = new TFile(file_down);
		auto down_graph = (TGraphAsymmErrors *)down_file->Get(name_file_down);

		for (int i = 0; i < down_graph->GetN(); i++)
		{
			double graph_xu = down_graph->GetErrorXhigh(i);
			double graph_xd = down_graph->GetErrorXlow(i);
			double graph_yu = down_graph->GetErrorYhigh(i);
			double graph_yd = down_graph->GetErrorYlow(i);
			double graph_x = down_graph->GetPointX(i);
			double graph_y = down_graph->GetPointY(i) * weight_file_down;
			down_graph->SetPoint(i, graph_x, graph_y);
			down_graph->SetPointError(i, graph_xd, graph_xu, graph_yd, graph_yu);
		}

		TFile *up_file = new TFile(file_up);
		auto up_graph = (TGraphAsymmErrors *)up_file->Get(name_file_up);

		for (int j = 0; j < up_graph->GetN(); j++)
		{
			double graph_xu = up_graph->GetErrorXhigh(j);
			double graph_xd = up_graph->GetErrorXlow(j);
			double graph_yu = up_graph->GetErrorYhigh(j);
			double graph_yd = up_graph->GetErrorYlow(j);
			double graph_x = up_graph->GetPointX(j);
			double graph_y = up_graph->GetPointY(j) * weight_file_up;
			up_graph->SetPoint(j, graph_x, graph_y);
			up_graph->SetPointError(j, graph_xd, graph_xu, graph_yd, graph_yu);
		}

		base_graph->GetListOfFunctions()->Clear();
		down_graph->GetListOfFunctions()->Clear();
		up_graph->GetListOfFunctions()->Clear();
		//////////////////////////////////

		// base_graph->SetLineWidth(2);
		base_graph->SetLineColor(kBlack);
		base_graph->SetMarkerColor(kBlack);
		base_graph->SetMarkerSize(2);
		base_graph->SetMarkerStyle(20);
		base_graph->GetYaxis()->SetTitleOffset(0.8);

		// down_graph->SetLineWidth(2);
		down_graph->SetLineColor(kRed);
		down_graph->SetMarkerColor(kRed);
		down_graph->SetMarkerSize(2);
		down_graph->SetMarkerStyle(20);

		// up_graph->SetLineWidth(2);
		up_graph->SetLineColor(kOrange);
		up_graph->SetMarkerColor(kOrange);
		up_graph->SetMarkerSize(2);
		up_graph->SetMarkerStyle(20);

		auto legend = new TLegend(0.54, 0.87, 0.87, 0.67);
		legend->AddEntry(base_graph, name_base, "lp");
		legend->AddEntry(down_graph, name_down, "lp");
		legend->AddEntry(up_graph, name_up, "lp");

		int nb_point = base_graph->GetN();
		TGraphAsymmErrors *syst_var_down = new TGraphAsymmErrors(nb_point);
		TGraphAsymmErrors *syst_var_up = new TGraphAsymmErrors(nb_point);

		for (int j = 0; j < nb_point; j++)
		{

			double base_xu = base_graph->GetErrorXhigh(j);
			double base_xd = base_graph->GetErrorXlow(j);
			if (base_xu == 0.0 && base_xd == 0.0)
			{
				base_xu = 0.1;
				base_xd = 0.1;
			}
			double base_x = base_graph->GetPointX(j);

			double base_y = base_graph->GetPointY(j);
			double down_y = down_graph->GetPointY(j);
			double up_y = up_graph->GetPointY(j);

			syst_var_down->SetPoint(j, base_x, 0.0);
			syst_var_down->SetPointError(j, base_xd, base_xu, (base_y - down_y) * 100. / (base_y), 0.0);

			syst_var_up->SetPoint(j, base_x, 0.0);
			syst_var_up->SetPointError(j, base_xd, base_xu, (base_y - up_y) * 100. / (base_y), 0.0);

			cout << base_y - up_y << endl;
		}

		if (combinaison_label == "int")
			updateGraphErrors(total_graph_int, syst_var_down, syst_var_up);
		if (combinaison_label == "bin1")
			updateGraphErrors(total_graph_bin1, syst_var_down, syst_var_up);
		if (combinaison_label == "bin2")
			updateGraphErrors(total_graph_bin2, syst_var_down, syst_var_up);
		if (combinaison_label == "bin3")
			updateGraphErrors(total_graph_bin3, syst_var_down, syst_var_up);

		syst_var_down->SetMaximum(max_y_syst);
		syst_var_down->SetMinimum(min_y_syst);
		syst_var_down->GetXaxis()->SetTitleSize(30);
		syst_var_down->GetXaxis()->SetTitleFont(43);
		syst_var_down->GetXaxis()->SetTitleOffset(4.);
		syst_var_down->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		syst_var_down->GetXaxis()->SetLabelSize(30);
		// Y axis ratio plot settings
		syst_var_down->SetTitle(Form(";%s;Variation (%)", var_label.Data()));
		syst_var_down->GetYaxis()->SetNdivisions(505);
		syst_var_down->GetYaxis()->SetTitleSize(30);
		syst_var_down->GetYaxis()->SetTitleFont(43);
		syst_var_down->GetYaxis()->SetTitleOffset(1.55);
		syst_var_down->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		syst_var_down->GetYaxis()->SetLabelSize(30);
		syst_var_down->SetLineColor(kRed);
		syst_var_down->SetMarkerColor(kRed);
		syst_var_down->SetFillColorAlpha(kRed, 0.5);
		// syst_var_down->SetFillStyle(4080);

		syst_var_up->SetLineColor(kOrange);
		syst_var_up->SetMarkerColor(kOrange);
		syst_var_up->SetFillColorAlpha(kOrange, 0.5);
		// syst_var_up->SetFillStyle(4080);

		TCanvas *cancG0 = new TCanvas("", "can0", 1500, 1000);
		cancG0->cd();
		float limit_lower_pad = 0.3;
		TPad *pad1 = new TPad("pad1", "pad1", 0, limit_lower_pad, 1, 1.0);
		pad1->SetBottomMargin(0.);
		// pad1->SetTicks(1, 1);

		pad1->Draw(); // Draw the upper pad: pad1
		pad1->cd();	  // pad1 becomes the current pad

		cout << log_label << endl;
		if (log_label == "log")
		{
			gPad->SetLogy();
			base_graph->SetMaximum(10.0 * base_graph->GetMaximum());
		}

		base_graph->Draw("AP ");
		down_graph->Draw("P");
		up_graph->Draw("P");

		legend->SetFillStyle(0);
		legend->SetLineWidth(0);
		legend->Draw("same ");

		// Labels
		double x_top_label = 0.90;

		TPaveText *CLAS12_Internal = new TPaveText(0.2, x_top_label, 0.288191, x_top_label + 0.1, "NDC");
		CLAS12_Internal->SetFillStyle(4050);
		CLAS12_Internal->SetLineColor(0);
		CLAS12_Internal->SetTextFont(42);
		CLAS12_Internal->SetTextSize(0.0599401);
		CLAS12_Internal->SetBorderSize(0);
		CLAS12_Internal->AddText(title_sys);
		CLAS12_Internal->Draw();

		TText *t = new TText();
		t->SetTextAlign(22);
		t->SetTextColorAlpha(kGray, 0.50);
		t->SetTextFont(40);
		t->SetTextSize(0.25);
		t->SetTextAngle(25);
		t->DrawTextNDC(.5, .53, "Preliminary");

		// lower plot will be in pad
		///////////////////////////////////////////////////////////////////////////////
		cancG0->cd(); // Go back to the main canvas before defining pad2
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, limit_lower_pad - 0.005);
		pad2->SetTopMargin(0);
		// pad2->SetTicks(1, 1);
		pad2->SetGridy();

		pad2->SetBottomMargin(0.3);
		pad2->Draw();
		pad2->cd();

		auto legend_1 = new TLegend(0.54, 0.87, 0.87, 0.67);
		// legend_1->AddEntry(ratio_hist, "Ratio V1/No Rad.", "lp");
		legend_1->SetFillStyle(0);
		legend_1->SetLineWidth(0);

		syst_var_down->Draw("a2");
		syst_var_up->Draw("e2");
		legend_1->Draw();
		pad2->Update();

		TString name_pdf = "Systematics_JPsi";
		if (i == 0)
			cancG0->SaveAs(name_pdf + ".pdf(");
		else if (i == (Systematics_array.size() - 1))
			cancG0->SaveAs(name_pdf + ".pdf)");
		else
			cancG0->SaveAs(name_pdf + ".pdf");
	}

	cout<<endl;
	cout<<"Integrated plot"<<endl;
	printGraphContent(total_graph_int);

	cout<<endl;
	cout<<"Bin 1"<<endl;
	printGraphContent(total_graph_bin1);

	cout<<endl;
	cout<<"Bin 2"<<endl;
	printGraphContent(total_graph_bin2);

	cout<<endl;
	cout<<"Bin 3"<<endl;
	printGraphContent(total_graph_bin3);

	gApplication->Terminate();

	return 0;
}
