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
#include <TApplication.h>
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

class BSA_Analysis : public Analysis
{
  public:

	//Setup histogrmms
	TH1D *nb_JPsi_Positive_helicity;
	TH1D *nb_JPsi_Negative_helicity;

	//class attributes for the analysis
	int bin_phi = 4;

	//New classes to handle the BSA binning and use only inbending 2018
	// Phi binning should be handled better in the future
    void Set_Binning_BSA_RGA()
	{
		// Store the binning and plotting configuration
		TString min_hist_MC = "2.0";
		TString min_hist = "2.6";
		TString max_hist = "3.5";
		TString bin_hist = "45";//"30.";//

		//labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0", "", "M2"});


		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==1 && phi>-180 && phi<-90", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==1 && phi>-90 && phi<0", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==1 && phi>0 && phi<90", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==1 && phi>90 && phi<180", "", "M2"});

		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==-1 && phi>-180 && phi<-90", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==-1 && phi>-90 && phi<0", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==-1 && phi>0 && phi<90", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && helicity==-1 && phi>90 && phi<180", "", "M2"});
		//shitty fast way to do the phi binning, need to rewrite this with a loop
		
	}

	void Set_Sample_to_RGA_inbending_2018()
	{
		Analysis_Sample.Setup_RGA_inbending_2018();
		RunGroup = "RGA";
	}
	
	void Setup_Histo_BSA()
	{
		nb_JPsi_Positive_helicity = new TH1D("nb_JPsi_Positive_helicity", "nb_JPsi_Positive_helicity", bin_phi, -180, 180);
		nb_JPsi_Positive_helicity->SetStats(kFALSE);
		nb_JPsi_Positive_helicity->GetYaxis()->SetTitleOffset(1.3);
		nb_JPsi_Positive_helicity->SetTitle(";Phi;nb JPsi");

		nb_JPsi_Negative_helicity = new TH1D("nb_JPsi_Negative_helicity", "nb_JPsi_Negative_helicity", bin_phi, -180, 180);
		nb_JPsi_Negative_helicity->SetStats(kFALSE);
		nb_JPsi_Negative_helicity->GetYaxis()->SetTitleOffset(1.3);
		nb_JPsi_Negative_helicity->SetTitle(";Phi;nb JPsi");
	}

	void Process_REC_BSA()
	{

			TCanvas *cancG0 = new TCanvas("", "can0", 1500, 1000);
			cancG0->cd();

		for (int i = 0; i < labels.size(); i++)
		{

			//////////////////////////////////////////////////
			// Set options for each label
			TString label = labels[i][0];
			TString xAxis_label = labels[i][1];
			TString min_histo_option = labels[i][2];
			TString max_histo_option = labels[i][3];
			TString nb_bins = labels[i][4];
			TString string_cut = labels[i][5];
			TString max_y = labels[i][6];
			TString output_string = labels[i][7];
			//////////////////////////////////////////////////

			cout << "\n";
			cout << "//////////////////////////////////////////////////"<< "\n";
			cout << "Doing " << i << " " << label << " plot"<< "\n";
			cout << "//////////////////////////////////////////////////"<< "\n";
			cout << "\n";

			TCut cut = string_cut.Data(); //""; //"1";////

			TString cut_string = cut.GetTitle();

			float maximum = (Analysis_Sample.filtered_Data_tree->GetMaximum(label));
			float minimum = (Analysis_Sample.filtered_Data_tree->GetMinimum(label));
			float min_histo_ini = stof((string)min_histo_option.Data());
			float max_histo_ini = stof((string)max_histo_option.Data());

			int nBins = stoi((string)nb_bins.Data());

			TH1D *Data_hist = new TH1D(Form("Data_hist_%i", i), "Data_hist", nBins, min_histo_ini, max_histo_ini);
			TH1D *Average_variable = new TH1D(Form("Average_variable_%i", i), "Average_variable", 100, 0., 11.);

			TCut weight_data = Form("%s", "weight");

			Analysis_Sample.filtered_Data_tree->Draw(label + ">>" + Form("Data_hist_%i", i), weight_data * data_cut * cut * exclusivity_cut * kinematic_cut);
			Analysis_Sample.filtered_Data_tree->Draw(variable + ">>" + Form("Average_variable_%i", i), weight_data * data_cut * cut * exclusivity_cut * kinematic_cut);
			average_variable.push_back(Average_variable->GetMean());
			sigma_variable.push_back(Average_variable->GetRMS());


			Data_hist->SetLineWidth(2);
			Data_hist->SetLineColor(kBlack);
			Data_hist->SetMarkerColor(kBlack);
			Data_hist->SetMarkerSize(2);
			Data_hist->SetMarkerStyle(20);
			Data_hist->SetTitle(";" + label + ";Events");
			Data_hist->SetStats(kFALSE);

			cout << "///////////////////////" << endl;
			cout << "FIT DATA" << endl;
			cout << "///////////////////////" << endl;
			Fit_Function Fit_func;
			Fit_func.Set_Data_hist(Data_hist);
			Fit_func.Set_Limits(min_fit, max_fit);
			int nb_fit = 0;

			while (!Fit_func.fit_status && nb_fit<40)
			{
				if (fit_procedure == "Default")
					Fit_func.Single_Gaussian_Int_fit("SLER", Form("func_%i", i));
				if (fit_procedure == "Crystall ball Pol 2 BG")
					Fit_func.Crystall_Ball_fit("SLER", Form("func_%i", i));
				if (fit_procedure == "Crystall ball exp BG")
					Fit_func.Crystall_Ball_fit_exp("SLR", Form("func_%i", i));
				if (fit_procedure == "Pol 2 BG")
					Fit_func.Single_Gaussian_Int_fit_Pol_BG_V2("SLER", Form("func_%i", i));
				if (fit_procedure == "Double Gaussian")
					Fit_func.Double_Gaussian_Fit("SELR", Form("func_%i", i));
				nb_fit++;
				cout << Fit_func.fit_status << endl;
			}

			double chi2 = Fit_func.chi2;
			double NDF = Fit_func.NDF;

			Fit_func.Draw_Functions();

			double nb_JPsi = Fit_func.Get_Integral_Signal();			 
			double error_nb_JPsi = Fit_func.Get_Integral_Error_Signal(); 
			nb_JPsi_Data.push_back(nb_JPsi);
			error_nb_JPsi_Data.push_back(error_nb_JPsi);
			cout << " nb_JPsi  " << nb_JPsi << endl;

			auto legend = new TLegend(0.54, 0.87, 0.90, 0.60);
			legend->AddEntry(Data_hist, Form("Data (%3.1f)", Data_hist->Integral()), "lp");
			legend->AddEntry(Fit_func.function_Signal, Form("J#psi fit (%3.1f #pm %3.1f)  ", nb_JPsi, error_nb_JPsi), "l");
			legend->AddEntry(Fit_func.function_Signal, Form("#Chi^{2}: %3.1f, NdF: %3.1f, #Chi^{2}/NdF: %3.1f ", chi2, NDF, chi2 / NDF), "");
			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");
			

			if (i == 0)
				cancG0->SaveAs(output_folder + name_pdf + ".pdf(");
			else
				cancG0->SaveAs(output_folder + name_pdf + ".pdf");
		}

	}

	void Finalize_BSA_Calculation()
	{

		for(int i=0; i<bin_phi; i++)
		{
			nb_JPsi_Positive_helicity->SetBinContent(i+1, nb_JPsi_Data[i]);
			nb_JPsi_Positive_helicity->SetBinError(i+1, error_nb_JPsi_Data[i]);
		}
		for(int i=0; i<bin_phi; i++)
		{
			nb_JPsi_Negative_helicity->SetBinContent(i+1, nb_JPsi_Data[i+bin_phi]);
			nb_JPsi_Negative_helicity->SetBinError(i+1, error_nb_JPsi_Data[i+bin_phi]);
		}

		TCanvas *cancG0 = new TCanvas("", "cancG0", 1500, 1000);
		cancG0->cd();
		nb_JPsi_Positive_helicity->Draw();
		cancG0->SaveAs(output_folder + name_pdf + ".pdf");

		TCanvas *cancG1 = new TCanvas("", "cancG1", 1500, 1000);
		cancG1->cd();
		nb_JPsi_Negative_helicity->Draw();
		cancG1->SaveAs(output_folder + name_pdf + ".pdf");

		TH1D *numerator = (TH1D *)nb_JPsi_Positive_helicity->Clone("nb_JPsi_Positive_helicity");
		numerator->Add(nb_JPsi_Negative_helicity, -1);
		TCanvas *cancG3 = new TCanvas("", "cancG3", 1500, 1000);
		cancG3->cd();
		numerator->Draw();
		cancG3->SaveAs(output_folder + name_pdf + ".pdf");

		TH1D *denominator = (TH1D *)nb_JPsi_Positive_helicity->Clone("nb_JPsi_Positive_helicity");
		denominator->Add(nb_JPsi_Negative_helicity);
		TCanvas *cancG4 = new TCanvas("", "cancG4", 1500, 1000);
		cancG4->cd();
		denominator->Draw();
		cancG4->SaveAs(output_folder + name_pdf + ".pdf");

		TH1D *BSA_hist = (TH1D *)numerator->Clone("nb_JPsi_Positive_helicity");
		BSA_hist->Divide(denominator);
		TCanvas *cancG2 = new TCanvas("", "cancG2", 1500, 1000);
		cancG2->cd();
		BSA_hist->SetTitle(";Phi;BSA");
		BSA_hist->Draw();
		cancG2->SaveAs(output_folder + name_pdf + ".pdf)");



	}

};

int BSA_Extraction()
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
	gRandom->SetSeed(0);

	cout<<"test"<<endl;

	///////////////////////////////////////////////////////
	////////////// Instantiate analysis class /////////////
	///////////////////////////////////////////////////////
	BSA_Analysis JPsi_BSA_analysis;

	cout<<"test"<<endl;

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
		JPsi_BSA_analysis.Set_Parameters(parameters_file);
	}
	else
	{
		parameters_file = argv[3];
		cout<<"Parameter file is: "<<parameters_file<<endl;
		JPsi_BSA_analysis.Set_Parameters(parameters_file);
	}

	// Allow to overide some parameters
	if (input.cmdOptionExists("-output_folder"))
	{
		JPsi_BSA_analysis.output_folder = input.getCmdOption("-output_folder");
	}

	///////////////////////////////////////////////////////
	////////////////   Setup the analysis  ////////////////
	///////////////////////////////////////////////////////
	//JPsi_BSA_analysis.Set_Sample_to_RGA();
	JPsi_BSA_analysis.Set_Sample_to_RGA_inbending_2018();
	JPsi_BSA_analysis.Set_Binning_BSA_RGA();

	JPsi_BSA_analysis.Setup_Histo_BSA();
	//JPsi_BSA_analysis.Setup_Latex_Table();
	///////////////////////////////////////////////////////

	////////////////////////////////////
	//////// Filter the samples ////////
	////////////////////////////////////
	JPsi_BSA_analysis.Analysis_Sample.Filter_Data_Tree(JPsi_BSA_analysis.kinematic_cut, JPsi_BSA_analysis.data_cut, JPsi_BSA_analysis.exclusivity_cut);
	time(&filtering_time);
	////////////////////////////////////

	/////////////////////////////////////////
	//////// Cross section calculation //////
	/////////////////////////////////////////
	JPsi_BSA_analysis.Process_REC_BSA();
	time(&fitting_time);
	//JPsi_BSA_analysis.Process_MC();
	//time(&MC_time);
	JPsi_BSA_analysis.Finalize_BSA_Calculation();

	//JPsi_BSA_analysis.Save_Histo_to_root();
	//JPsi_BSA_analysis.Latex_Table.Format();
	//JPsi_BSA_analysis.Latex_Table.Save();
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

	cout << " -- Terminated -- " << endl;

	//gDebug = 2;

	// gROOT->GetListOfHistograms()->Delete();
	gApplication->Terminate();

	// exit(0);
	// gSystem->Exit(0);
	// gSystem->ProcessEvents();  // Ensure pending events are processed
	// cout << " Here after process event " << endl;
	// cout << " .q line " << endl;
	// gROOT->ProcessLine(".q");
	// gSystem->Exit(0);          // Forceful termination if necessary

	return 0;
}
