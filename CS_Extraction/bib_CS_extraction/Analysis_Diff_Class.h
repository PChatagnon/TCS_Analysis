#ifndef Analysis_Diff_Class
#define Analysis_Diff_Class

#include "Utils.h"
#include "Analysis_Class.h"

class Analysis_Diff : public Analysis
{
public:
	double Eg_min = 8.2;  // //8.205; //
	double Eg_max = 9.28; //  // 10.6;  //

	std::vector<double> average_variable_2{};
	std::vector<double> sigma_variable_2{};

	TString variable_2 = "Epho";

	////// Store histogramms and Graphs //////
	TGraphAsymmErrors JPsi_CS_Graph_C;
	/////////////////////////////////

	Analysis_Diff()
	{
		bin_id = "bin 1";
		variable = "-t";
		output_folder = "/Users/pc281174/Desktop/JPsi_analysis/CS_Extraction/t_cross_section/";
		kinematic_cut = "pass_EC_cut &&  Proton.Theta()*180./3.141592<35. && M>2.7 && (Electron.P() > 1.7) && (Positron.P() > 1.7) && positron_SF>0.15 && electron_SF>0.15 && ( Positron.P()<4.0 || (Positron.P()>4.0 && positron_score>0.05)) && ( Electron.P()<4.0 || (Electron.P()>4.0 && electron_score>0.05))  && positron_HTCC_ECAL_match==1. && electron_HTCC_ECAL_match==1.";
		kinematic_cut_BG = "Proton.Theta()*180./3.141592<35. &&  M>2.7 && (Electron.P() > 1.7) && (Positron.P() > 1.7)";
	}

	//////  Setup binnings  //////
	void Set_Binning_t_diff_1()
	{

		Eg_min = 8.2;
		Eg_max = 9.28;

		bin_id = "bin 1";

		TString min_hist_MC = "1.8";
		TString min_hist = "2.7";
		TString max_hist = "3.3";
		TString bin_hist = "40.";

		labels = {};
		labels_MC = {};

		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>0.77 && -t<1.00 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.00 && -t<1.5 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.5 && -t<2.0 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>2.0 && -t<2.5 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>2.5 && -t<4.5 ", Eg_min, Eg_max), "", "M2"});

		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>0.77 && -t_Gen<1.00 ", Eg_min, Eg_max), "50", "0.77", "1.00"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.00 && -t_Gen<1.5 ", Eg_min, Eg_max), "50", "1.00", "1.5"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.5 && -t_Gen<2.0 ", Eg_min, Eg_max), "50", "1.5", "2.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>2.0 && -t_Gen<2.5 ", Eg_min, Eg_max), "50", "2.0", "2.5"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>2.5 && -t_Gen<4.5 ", Eg_min, Eg_max), "50", "2.5", "4.5"});

		name_pdf = Form(name_pdf + "%.2f_%.2f", Eg_min, Eg_max);
	}

	void Set_Binning_t_diff_2()
	{

		Eg_min = 9.28;
		Eg_max = 10.0;

		bin_id = "bin 2";

		TString min_hist_MC = "1.8";
		TString min_hist = "2.7";
		TString max_hist = "3.3";
		TString bin_hist = "40.";

		labels = {};
		labels_MC = {};

		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>0.5 && -t<0.75 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>0.75 && -t<1.0 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.0 && -t<1.25 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.25 && -t<1.5 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.5 && -t<1.75 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.75 && -t<2.0 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>2.0 && -t<2.5 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>2.5 && -t<3.0 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>3.0 && -t<6.0 ", Eg_min, Eg_max), "", "M2"});

		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>0.5 && -t_Gen<0.75 ", Eg_min, Eg_max), "50", "0.5", "0.75"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>0.75 && -t_Gen<1.0 ", Eg_min, Eg_max), "50", "0.75", "1.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.0 && -t_Gen<1.25 ", Eg_min, Eg_max), "50", "1.0", "1.25"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.25 && -t_Gen<1.5 ", Eg_min, Eg_max), "50", "1.25", "1.5"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.5 && -t_Gen<1.75 ", Eg_min, Eg_max), "50", "1.5", "1.75"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.75 && -t_Gen<2.0 ", Eg_min, Eg_max), "50", "1.75", "2.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>2.0 && -t_Gen<2.5 ", Eg_min, Eg_max), "50", "2.0", "2.5"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>2.5 && -t_Gen<3.0 ", Eg_min, Eg_max), "50", "2.5", "3.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>3.0 && -t_Gen<6.0 ", Eg_min, Eg_max), "50", "3.0", "6.0"});

		name_pdf = Form(name_pdf + "%.2f_%.2f", Eg_min, Eg_max);
	}

	void Set_Binning_t_diff_3()
	{

		Eg_min = 10.0;
		Eg_max = 10.6;

		bin_id = "bin 2";

		TString min_hist_MC = "1.8";
		TString min_hist = "2.7";
		TString max_hist = "3.3";
		TString bin_hist = "40.";

		labels = {};
		labels_MC = {};

		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>0.5 && -t<0.75 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>0.75 && -t<1.0 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.0 && -t<1.25 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>1.25 && -t<2.0 ", Eg_min, Eg_max), "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, Form("status_prot<4000 && M>2.0 && Epho>%f && Epho<%f && -t>2.0 && -t<6.0 ", Eg_min, Eg_max), "", "M2"});

		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>0.5 && -t_Gen<0.75 ", Eg_min, Eg_max), "50", "0.5", "0.75"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>0.75 && -t_Gen<1.0 ", Eg_min, Eg_max), "50", "0.75", "1.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.0 && -t_Gen<1.25 ", Eg_min, Eg_max), "50", "1.0", "1.25"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>1.25 && -t_Gen<2.0 ", Eg_min, Eg_max), "50", "1.25", "2.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, Form("weight<100 && M_Gen_2>2.0 && Epho_Gen>%f && Epho_Gen<%f && -t_Gen>2.0 && -t_Gen<6.0 ", Eg_min, Eg_max), "50", "2.0", "6.0"});
		
		name_pdf = Form(name_pdf + "%.2f_%.2f", Eg_min, Eg_max);
	}
	//////////////////////////////

	////// Method for Latex Tables //////
	void Setup_Latex_Table()
	{
		Latex_Table = Latex_Table_writter("", output_folder.Data(), "t");
		Latex_Table.Set_output_name(((string)name_pdf.Data()) + "_Latex_Table.txt");
		Latex_Table.mode = "2D";
	}
	//////////////////////////////

	/////// Processing function for the integrated CS //////
	void Process_REC_Diff()
	{

		JPsi_CS_Graph.Set(labels.size());

		// TString cut_Eg_str =
		TCut cut_Eg = Form("status_prot<4000 && Epho>%f && Epho<%f ", Eg_min, Eg_max);
		TH1D *Average_Eg = new TH1D("Average_Eg", "Average_Eg", 100, 0., 11.);
		Analysis_Sample.filtered_Data_tree->Draw("Epho>>Average_Eg", cut_Eg);
		cout << "Average Eg " << Average_Eg->GetMean() << " sigma " << Average_Eg->GetRMS() << endl;

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
			cout << "//////////////////////////////////////////////////"
				 << "\n";
			cout << "Doing " << label << " plot"
				 << "\n";
			cout << "//////////////////////////////////////////////////"
				 << "\n";
			cout << "\n";

			TCut cut = string_cut.Data(); //""; //"1";////

			TString cut_string = cut.GetTitle();

			float maximum = (Analysis_Sample.filtered_Data_tree->GetMaximum(label));
			float minimum = (Analysis_Sample.filtered_Data_tree->GetMinimum(label));

			float min_histo_ini = stof((string)min_histo_option.Data());
			float max_histo_ini = stof((string)max_histo_option.Data());

			int nBins = stoi((string)nb_bins.Data());

			TH1D *Data_hist = new TH1D("Data_hist", "Data_hist", nBins, min_histo_ini, max_histo_ini);
			TH2D *Average_variable = new TH2D("Average_variable", "Average_variable", 100, 0., 11. , 100, 0., 11.);

			// TCanvas *cancG1 = new TCanvas("", "can1",1500,1000);
			TCut weight_data = Form("%s", "weight");
			// TCut weight_data = Form("%s", "1");

			Analysis_Sample.filtered_Data_tree->Draw(label + ">>Data_hist", weight_data * data_cut * cut * exclusivity_cut * kinematic_cut);
			Analysis_Sample.filtered_Data_tree->Draw(variable + ":" + variable_2 + + ">>Average_variable", weight_data * data_cut * cut * exclusivity_cut * kinematic_cut);
			average_variable.push_back(Average_variable->ProjectionY()->GetMean());
			sigma_variable.push_back(Average_variable->ProjectionY()->GetRMS());

			average_variable_2.push_back(Average_variable->ProjectionX()->GetMean());
			sigma_variable_2.push_back(Average_variable->ProjectionX()->GetRMS());

			cout<<"Mean and RMS of variables "<<Average_variable->ProjectionX()->GetMean()<<" "<<Average_variable->ProjectionY()->GetMean()<<endl;

			Data_hist->SaveAs("plot.pdf");

			int nbBins = Data_hist->GetNbinsX();
			float min_histo = Data_hist->GetXaxis()->GetXmin();
			float max_histo = Data_hist->GetXaxis()->GetXmax();
			if (debug)
				cout << "Bining and range " << nbBins << " " << min_histo << " " << max_histo << "\n";

			double r = 0.0;
			r = Data_hist->GetBinContent(1) + Data_hist->GetBinContent(0);
			Data_hist->SetBinContent(1, r);
			r = Data_hist->GetBinContent(nBins) + Data_hist->GetBinContent(nBins + 1);
			Data_hist->SetBinContent(nBins, r);

			Data_hist->SetLineWidth(2);
			Data_hist->SetLineColor(kBlack);
			Data_hist->SetMarkerColor(kBlack);
			Data_hist->SetMarkerSize(2);
			Data_hist->SetMarkerStyle(20);
			Data_hist->SetTitle(";" + label + ";Events");
			Data_hist->SetStats(kFALSE);

			auto legend = new TLegend(0.54, 0.87, 0.90, 0.60);
			legend->AddEntry(Data_hist, Form("Data (%3.1f)", Data_hist->Integral()), "lp");

			THStack *hs = new THStack("hs", "");
			THStack *hs_JPsi = new THStack("hs_JPsi", "");

			TH1D *BG_hist = new TH1D("BG_hist", "BG_hist", nBins, min_histo, max_histo);

			// TTree *small = ntuple->CopyTree("py>2");

			for (int j = 0; j < Analysis_Sample.samples.size(); j++)
			{

				if (Analysis_Sample.samples[j][3] == "Jpsi_Rad_Corr")
					continue;

				int nbEvents_sample = Analysis_Sample.ngen[j];
				if (debug)
					cout << "Nb of events in " << Analysis_Sample.samples[j][0] << " : " << nbEvents_sample << "\n";

				// Xsec BG
				double xsec = stof((string)Analysis_Sample.samples[j][1].Data());
				double lumi_sample = stof((string)Analysis_Sample.samples[j][4].Data());

				TString hist_name = Form("sample_hist_%s_%i_%i", Analysis_Sample.samples[j][3].Data(), j, i);
				TH1D *sample_hist = new TH1D(hist_name, "", nbBins, min_histo, max_histo);
				TCut weight = Form("%s*%f*%f*%f/(%i)", "weight", xsec, Analysis_Sample.lumi_factor, lumi_sample, nbEvents_sample);

				// if (Analysis_Sample.samples[j][3] != "BH TCSGen")weight = Form("%s*%s*%f*%f*%f/(%i*%s)", "weight", "virtual_flux_Frixione_Gen", xsec, Analysis_Sample.lumi_factor, lumi_sample, nbEvents_sample, "virtual_flux_Gen");

				// TString hist_Acc = Form("sample_Acc_%s", Analysis_Sample.samples[j][3].Data());
				// TH1F *sample_Acc = new TH1F(hist_Acc, "", nbBins, min_histo, max_histo);
				// TCut weight_Acc = Form("%s*%f", "weight", xsec);
				//  TCut weight_Acc = Form("%s*%s*%f/(%s)", "weight", "virtual_flux_Frixione_Gen", xsec, "virtual_flux_Gen");

				Analysis_Sample.reduced_samples_tree[j]->Draw(label + ">>" + hist_name, cut * weight);
				// reduced_samples_tree[j]->Draw(label + ">>" + hist_Acc, cut * weight_Acc);

				if (debug)
					cout << cut * weight << "\n";
				if (debug)
					cout << "Number of entries " << sample_hist->GetEntries() << "\n";

				// UnderFlow-OverFlow
				r = sample_hist->GetBinContent(1) + sample_hist->GetBinContent(0);
				sample_hist->SetBinContent(1, r);
				r = sample_hist->GetBinContent(nBins) + sample_hist->GetBinContent(nBins + 1);
				sample_hist->SetBinContent(nBins, r);

				sample_hist->SetLineWidth(0);
				sample_hist->SetLineColor(kBlack);
				sample_hist->SetMarkerSize(0);
				sample_hist->SetFillColorAlpha(std::stoi(Analysis_Sample.samples[j][2].Data()), 0.65);
				sample_hist->SetStats(kFALSE);

				hs->SetTitle(";" + label + ";Events");
				hs->Add(sample_hist);

				hs_JPsi->SetTitle(";" + label + ";Events");
				if (Analysis_Sample.samples[j][3] == "J#psi")
					hs_JPsi->Add(sample_hist);

				if (debug)
					cout << Analysis_Sample.samples[j][0] << "\n";
				BG_hist->Add(sample_hist);

				// legend->AddEntry(sample_hist, Form("%s (%3.1f)", Analysis_Sample.samples[j][3].Data(), sample_hist->Integral()), "f1");
			}

			// Ratio plots and uncertainty
			TH1D *ratio_hist = (TH1D *)Data_hist->Clone("Data_hist");
			ratio_hist->Divide(BG_hist);

			TH1D *ratio_hist_uncertainty = (TH1D *)Data_hist->Clone("Data_hist");
			for (int ii = 1; ii < nBins + 1; ii++)
			{
				ratio_hist_uncertainty->SetBinContent(ii, 1);
				if (BG_hist->GetBinContent(ii) > 0)
				{
					ratio_hist_uncertainty->SetBinError(ii, BG_hist->GetBinError(ii) / BG_hist->GetBinContent(ii));
				}
				else
				{
					ratio_hist_uncertainty->SetBinError(ii, 1);
				}
			}

			TCanvas *cancG0 = new TCanvas("", "can0", 1500, 1000);
			cancG0->cd();

			float limit_lower_pad = 0.0;
			if (ratio_pad || pull_pad)
				limit_lower_pad = 0.3;
			// Upper plot will be in pad1
			TPad *pad1 = new TPad("pad1", "pad1", 0, limit_lower_pad, 1, 1.0);

			if (ratio_pad || pull_pad)
				pad1->SetBottomMargin(0.003);

			float max_display = (hs->GetMaximum()) * 1.5;
			hs->SetMaximum(max_display);
			// hs_JPsi->SetMaximum(max_display);

			pad1->Draw(); // Draw the upper pad: pad1
			pad1->cd();	  // pad1 becomes the current pad

			if (debug)
				cout << " max _ y " << max_y << "\n";
			if (max_y != "")
			{
				float max_y_histo = std::stof(max_y.Data());

				Data_hist->SetMaximum(max_y_histo);
			}
			else
				Data_hist->SetMaximum(Data_hist->GetMaximum() * 1.6);

			Data_hist->Draw("e");

			Fit_Function Fit_func;
			Fit_func.Set_Data_hist(Data_hist);
			Fit_func.Set_Limits(min_fit, max_fit);
			
			if(fit_procedure=="Default")
				Fit_func.Single_Gaussian_Int_fit("SLER", Form("func_%i", i));

			if(fit_procedure=="Crystall ball Pol 2 BG")
				Fit_func.Crystall_Ball_fit("SLER", Form("func_%i", i));

			if(fit_procedure=="Crystall ball exp BG")
				Fit_func.Crystall_Ball_fit_exp("SLER", Form("func_%i", i));

			if(fit_procedure=="Pol 2 BG")
				Fit_func.Single_Gaussian_Int_fit_Pol_BG_V2("SLER", Form("func_%i", i));

			if(fit_procedure=="Double Gaussian")
				Fit_func.Double_Gaussian_Fit("SLR",Form("func_%i", i));

			double chi2 = Fit_func.chi2;
			double NDF = Fit_func.NDF;

			if (ratio_pad)
			{
				hs->Draw("e hist same");
				BG_hist->SetMarkerSize(0);
				BG_hist->SetFillColor(kGray);
				BG_hist->SetFillStyle(3144);
				BG_hist->Draw("same  e2");
				Data_hist->Draw("e same");
			}

			Fit_func.Draw_Functions();

			double nb_JPsi = Fit_func.Get_Integral_Signal();			 //(JPsi_signal->Integral(0., 10.)) / (Data_hist->GetXaxis()->GetBinWidth(2));
			double error_nb_JPsi = Fit_func.Get_Integral_Error_Signal(); // sqrt(nb_JPsi); // sigma_integral;//nb_JPsi * sqrt((error_amp_fit * error_amp_fit) / (amp_fit * amp_fit) + (error_sigma_fit * error_sigma_fit) / (sigma_fit * sigma_fit) + 2. * covMatrix(0, 2) / (sigma_fit * amp_fit));

			cout << " nb_JPsi  " << nb_JPsi << endl;

			nb_JPsi_Data.push_back(nb_JPsi);
			error_nb_JPsi_Data.push_back(error_nb_JPsi);

			double nb_JPsi_C = (Data_hist->Integral(Data_hist->FindBin(3.0), Data_hist->FindBin(10.)));
			nb_JPsi_Data_C.push_back(nb_JPsi_C);

			w_c.push_back(w_c_from_BG_estimation); // Data_hist->Integral(Data_hist->FindBin(2.5), Data_hist->FindBin(2.9)) / BG_hist->Integral(BG_hist->FindBin(2.5), BG_hist->FindBin(2.9)));

			legend->AddEntry(Fit_func.function_Signal, Form("J#psi fit (%3.1f #pm %3.1f)  ", nb_JPsi, error_nb_JPsi), "l");
			legend->AddEntry(Fit_func.function_Signal, Form("#Chi^{2} %3.1f, NdF %3.1f, #Chi^{2}/NdF %3.1f ", chi2, NDF, chi2 / NDF), "");
			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");

			TText *t = new TText();
			t->SetTextAlign(22);
			t->SetTextColorAlpha(kGray, 0.50);
			t->SetTextFont(40);
			t->SetTextSize(0.25);
			t->SetTextAngle(25);

			// Labels
			double x_top_label = 0.90;

			TPaveText *CLAS12_Internal = new TPaveText(0.10, x_top_label, 0.288191, x_top_label + 0.1, "NDC");
			CLAS12_Internal->SetFillStyle(4050);
			CLAS12_Internal->SetLineColor(0);
			CLAS12_Internal->SetTextFont(42);
			CLAS12_Internal->SetTextSize(0.0599401);
			CLAS12_Internal->SetBorderSize(0);
			CLAS12_Internal->AddText(Form("CLAS12 Preliminary"));
			CLAS12_Internal->Draw();

			t->DrawTextNDC(.5, .53, "Preliminary");

			///////////////////////////////////////////////////////////////////////////////
			// lower plot will be in pad
			///////////////////////////////////////////////////////////////////////////////
			if (ratio_pad)
			{
				cancG0->cd(); // Go back to the main canvas before defining pad2
				TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, limit_lower_pad);
				pad2->SetTopMargin(0);
				// pad2->SetTicks(1, 1);
				pad2->SetGridy();

				pad2->SetBottomMargin(0.3);
				// pad2->SetGridx(); // vertical grid
				pad2->Draw();
				pad2->cd(); // pad2 becomes the current pad

				ratio_hist_uncertainty->SetFillColor(42);
				// ratio_hist_uncertainty->SetFillStyle(3001);
				ratio_hist_uncertainty->SetLineColor(1);
				ratio_hist_uncertainty->SetLineWidth(1);
				ratio_hist_uncertainty->SetMarkerSize(0);

				ratio_hist_uncertainty->SetMaximum(2.0);
				ratio_hist_uncertainty->SetMinimum(0.0);
				ratio_hist_uncertainty->Draw("e2");
				ratio_hist->Draw("ep same");
				pad2->Update();

				// Ratio plot (h3) settings
				ratio_hist_uncertainty->SetTitle(""); // Remove the ratio title

				// Y axis ratio plot settings
				ratio_hist_uncertainty->GetYaxis()->SetTitle("Data/MC ratio");
				ratio_hist_uncertainty->GetXaxis()->SetTitle(xAxis_label);
				ratio_hist_uncertainty->GetYaxis()->SetNdivisions(505);
				ratio_hist_uncertainty->GetYaxis()->SetTitleSize(30);
				ratio_hist_uncertainty->GetYaxis()->SetTitleFont(43);
				ratio_hist_uncertainty->GetYaxis()->SetTitleOffset(1.55);
				ratio_hist_uncertainty->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
				ratio_hist_uncertainty->GetYaxis()->SetLabelSize(30);

				// X axis ratio plot settings
				ratio_hist_uncertainty->GetXaxis()->SetTitleSize(30);
				ratio_hist_uncertainty->GetXaxis()->SetTitleFont(43);
				ratio_hist_uncertainty->GetXaxis()->SetTitleOffset(4.);
				ratio_hist_uncertainty->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
				ratio_hist_uncertainty->GetXaxis()->SetLabelSize(30);

				auto legendR = new TLegend(0.81, 0.95, 0.9, 0.75);
				legendR->AddEntry(ratio_hist_uncertainty, "MC Uncert.", "f1");
				legendR->Draw("same ");
				legendR->SetFillStyle(0);
				legendR->SetLineWidth(0);
			}

			if (pull_pad)
			{
				cancG0->cd(); // Go back to the main canvas before defining pad2
				TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, limit_lower_pad);
				pad2->SetTopMargin(0);
				pad2->SetGridy();
				pad2->SetBottomMargin(0.3);
				pad2->Draw();
				pad2->cd();


				TH1F *pulls_histo = new TH1F("pullsHistogram", "", Data_hist->GetNbinsX(), Data_hist->GetXaxis()->GetXmin(), Data_hist->GetXaxis()->GetXmax());
				for (int i = 1; i <= Data_hist->GetNbinsX(); ++i)
				{
					double dataValue = Data_hist->GetBinContent(i);
					double dataError = Data_hist->GetBinError(i);
					double fitValue = Fit_func.function->Eval(Data_hist->GetBinCenter(i));
					double pull = (dataValue - fitValue) / dataError;
					if(dataError==0.0)
						pull = 0.0;
					//cout<<"pull "<<pull<<endl;
					pulls_histo->SetBinContent(i,pull);
				}

				pulls_histo->SetStats(kFALSE);
				pulls_histo->SetFillColor(42);
				pulls_histo->SetLineColor(1);
				pulls_histo->SetLineWidth(0);
				pulls_histo->SetMarkerSize(0);
				pulls_histo->SetMaximum(5.0);
				pulls_histo->SetMinimum(-5.0);
				pulls_histo->Draw("hist");
				pulls_histo->SaveAs("pull.root");
				pad2->Update();

				pulls_histo->SetTitle("");
				pulls_histo->GetYaxis()->SetTitle("Pulls");
				pulls_histo->GetXaxis()->SetTitle(xAxis_label);
				pulls_histo->GetYaxis()->SetNdivisions(505);
				pulls_histo->GetYaxis()->SetTitleSize(30);
				pulls_histo->GetYaxis()->SetTitleFont(43);
				pulls_histo->GetYaxis()->SetTitleOffset(1.55);
				pulls_histo->GetYaxis()->SetLabelFont(43);
				pulls_histo->GetYaxis()->SetLabelSize(30);

				pulls_histo->GetXaxis()->SetTitleSize(30);
				pulls_histo->GetXaxis()->SetTitleFont(43);
				pulls_histo->GetXaxis()->SetTitleOffset(4.);
				pulls_histo->GetXaxis()->SetLabelFont(43);
				pulls_histo->GetXaxis()->SetLabelSize(30);
			}
			///////////////////////////////////////////////////////////////////////////////

			if (i == 0)
				cancG0->SaveAs(output_folder + name_pdf + ".pdf(");
			else
				cancG0->SaveAs(output_folder + name_pdf + ".pdf");
			cancG0->SaveAs(output_string + ".pdf");
			cancG0->SaveAs(output_string + ".png");

			

			/////////////////////////////////////////
			/////Acc calculation
			/////////////////////////////////////////

			/////Simulate the background and normalize simulations to data
			float nb_event_bg = (Fit_func.function_BG->Integral(min_histo, max_histo)) / (Data_hist->GetXaxis()->GetBinWidth(2));
			TString name_BG_func = Form("func_%i_bg_func", i);

			TH1D *h_only_JPsi = (TH1D *)(hs_JPsi->GetStack()->Last());
			double nb_JPsi_integral_MC_raw = h_only_JPsi->Integral(h_only_JPsi->FindBin(2.95), h_only_JPsi->FindBin(3.15));

			double additional_norm_factor = 4.0;
			double normalization_MC_to_data = additional_norm_factor * nb_JPsi; // nb_JPsi_integral_MC_raw;
			cout << "Normalization factor for signal in MC: " << normalization_MC_to_data << endl;

			int iteration_Acc = 1000;
			double nb_jpsi_per_iteration[iteration_Acc];
			Fit_Function Fit_func_MC;

			TCanvas *cancAcc = new TCanvas("", "can0", 1500, 1000);

			TH1D *BG_add_hist = new TH1D(Form("BG_add_hist_%i", i), "", nbBins, min_histo, max_histo);
			TH1D *h_JPsi_norm = new TH1D(Form("h_JPsi_norm_%i", i), "", nbBins, min_histo, max_histo);
			TH1D *hlast = new TH1D(Form("h_last_%i", i), "", nbBins, min_histo, max_histo);
			//TH1D *hs_JPsi_BG = new TH1D(Form("h_last_%i", i), "", nbBins, min_histo, max_histo);
			THStack *hs_JPsi_BG = new THStack("hs_JPsi_BG", "");

			for (int l = 0; l < iteration_Acc; l++)
			{
				BG_add_hist->Reset();
				h_JPsi_norm->Reset();
				hlast->Reset();

				BG_add_hist->FillRandom(name_BG_func, nb_event_bg);
				BG_add_hist->SetLineWidth(0);
				BG_add_hist->SetLineColor(kBlack);
				BG_add_hist->SetMarkerSize(0);
				BG_add_hist->SetFillColorAlpha(kCyan, 0.65);
				BG_add_hist->SetStats(kFALSE);

				h_JPsi_norm->SetFillColorAlpha(420, 0.50);
				h_JPsi_norm->SetLineWidth(0);
				for (int n = 0; n < normalization_MC_to_data; n++)
				{
					h_JPsi_norm->Fill(h_only_JPsi->GetRandom());
				}

				// h_only_JPsi->Scale(normalization_MC_to_data);
				double nb_JPsi_integral_MC = h_JPsi_norm->Integral();

				
				/////End simulate the background

				hlast->Add(h_JPsi_norm);
				hlast->Add(BG_add_hist);
				


				// cout << "///////////////////////" << endl;
				// cout << "FIT MC" << endl;
				/// cout << "///////////////////////" << endl;
				
				redefineErrors(hlast);
				Fit_func_MC.Set_Data_hist(hlast);
				Fit_func_MC.Set_Limits(min_fit, max_fit);
				// Fit_func_MC.Single_Gaussian_fit("SLER", Form("func_MC_%i", i));

				if (fit_procedure == "Default")
					Fit_func_MC.Single_Gaussian_Int_fit("QSLER", Form("func_MC_%i", i));

				if (fit_procedure == "Crystall ball Pol 2 BG")
					Fit_func_MC.Crystall_Ball_fit("QSLER", Form("func_MC_%i", i));

				if (fit_procedure == "Crystall ball exp BG")
					Fit_func_MC.Crystall_Ball_fit_exp("QSLER", Form("func_MC_%i", i));

				if (fit_procedure == "Pol 2 BG")
					Fit_func_MC.Single_Gaussian_Int_fit_Pol_BG_V2("QSLER", Form("func_MC_%i", i));

				if (fit_procedure == "Double Gaussian")
					Fit_func_MC.Double_Gaussian_Fit("QSLR", Form("func_MC_%i", i));

				double nb_JPsi_MC = Fit_func_MC.Get_Integral_Signal(); // nb_JPsi_integral; //
				nb_jpsi_per_iteration[l] = nb_JPsi_MC;

				if (l == 0)
				{
					cancAcc->cd();

					hlast->SetTitle(";" + label + ";Events");
					hlast->GetYaxis()->SetTitleOffset(1.3);
					hlast->SetStats(kFALSE);
					hlast->Draw("");

					hlast->SetLineWidth(1);
					hlast->SetLineColor(kBlack);
					hlast->SetMarkerStyle(1);
					hlast->SetMarkerSize(1);

					hs_JPsi_BG->Add(h_JPsi_norm);
					hs_JPsi_BG->Add(BG_add_hist);

					hlast->Draw("");
					hs_JPsi_BG->Draw("e same hist");
					Fit_func_MC.Draw_Functions();

					auto legend_acc = new TLegend(0.54, 0.87, 0.90, 0.60);
					legend_acc->AddEntry(h_only_JPsi, Form("nb JPsi %3.1f", nb_JPsi_integral_MC), "f1");
					legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("J#psi fit (%3.1f #pm %3.1f)  ", Fit_func_MC.Get_Integral_Signal(), Fit_func_MC.Get_Integral_Error_Signal()), "l");
					legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("#Chi^{2} %3.1f, NdF %3.1f, #Chi^{2}/NdF %3.1f ", Fit_func_MC.chi2, Fit_func_MC.NDF, Fit_func_MC.chi2 / Fit_func_MC.NDF), "");
					legend_acc->SetFillStyle(0);
					legend_acc->SetLineWidth(0);
					legend_acc->Draw("same ");

					cancAcc->SaveAs(output_folder + name_pdf + ".pdf");
				}
			}

			double numerator_Acc = 0;
			for (int i = 0; i < iteration_Acc; ++i)
			{
				numerator_Acc += nb_jpsi_per_iteration[i];
			}

			double iteration_in_double = static_cast<double>(iteration_Acc);
			cout << "total rec " << numerator_Acc << endl;
			cout << "iteration_in_double " << iteration_in_double << endl;
			cout << "Acc num " << nb_JPsi_integral_MC_raw * numerator_Acc / (additional_norm_factor * nb_JPsi * iteration_in_double) << endl;
			Acc_Num.push_back(nb_JPsi_integral_MC_raw * numerator_Acc / (additional_norm_factor * nb_JPsi * iteration_in_double)); //(sample_Acc->Integral());

			TCanvas *cancAccIteration = new TCanvas("", "cancAccIteration", 1500, 1000);
			cancAccIteration->cd();
			TH1D *histAccIteration = new TH1D(Form("histAccIteration%i", i), "", 50, 0.0 , 2.0 * numerator_Acc / iteration_in_double);
			for (int i = 0; i < iteration_Acc; ++i)
			{
				histAccIteration->Fill(nb_jpsi_per_iteration[i]);
			}

			gStyle->SetOptStat(1110);
			histAccIteration->Draw("hist");
			histAccIteration->SetTitle(";Fitted events;Counts");
			
			cancAccIteration->SaveAs(output_folder + name_pdf + ".pdf");

			gStyle->SetOptStat(1);


			/*
			/////////////////////////////////////////
			/////Acc calculation
			/////////////////////////////////////////

			/////Simulate the background
			float nb_event_bg = (Fit_func.function_BG->Integral(min_histo, max_histo)) / (Data_hist->GetXaxis()->GetBinWidth(2));
			TH1D *BG_add_hist = new TH1D("BG_add_hist", "", nbBins, min_histo, max_histo);
			TString name_BG_func = Form("func_%i_bg_func", i);
			BG_add_hist->FillRandom(name_BG_func, 50000);
			BG_add_hist->Scale(nb_event_bg / 50000.);
			// BG_add_hist->Scale(nb_event_bg/50000.); //Remove the BH background
			BG_add_hist->SetLineWidth(0);
			BG_add_hist->SetLineColor(kBlack);
			BG_add_hist->SetMarkerSize(0);
			BG_add_hist->SetFillColorAlpha(kCyan, 0.65);
			BG_add_hist->SetStats(kFALSE);
			TH1D *h_only_JPsi = (TH1D *)(hs_JPsi->GetStack()->Last());
			double nb_JPsi_integral = h_only_JPsi->Integral();
			hs_JPsi->Add(BG_add_hist);
			/////End simulate the background

			TCanvas *cancAcc = new TCanvas("", "can0", 1500, 1000);
			cancAcc->cd();
			TH1D *hlast = (TH1D *)(hs_JPsi->GetStack()->Last());
			hlast->SetTitle(";" + label + ";Events");
			hlast->GetYaxis()->SetTitleOffset(1.3);
			hlast->Draw("hist");
			hs_JPsi->Draw("e hist same");

			Fit_Function Fit_func_MC;
			Fit_func_MC.Set_Data_hist(hlast);
			Fit_func_MC.Set_Limits(min_fit, max_fit);

			if(fit_procedure=="Default")
				Fit_func_MC.Single_Gaussian_Int_fit("SLER", Form("func_MC_%i", i));

			if(fit_procedure=="Crystall ball Pol 2 BG")
				Fit_func_MC.Crystall_Ball_fit("SLER", Form("func_MC_%i", i));

			if(fit_procedure=="Crystall ball exp BG")
				Fit_func_MC.Crystall_Ball_fit_exp("SLER", Form("func_MC_%i", i));

			if(fit_procedure=="Pol 2 BG")
				Fit_func_MC.Single_Gaussian_Int_fit_Pol_BG_V2("SLER", Form("func_MC_%i", i));

			if(fit_procedure=="Double Gaussian")
				Fit_func_MC.Double_Gaussian_Fit("SLR",Form("func_MC_%i", i));

			hlast->Draw("hist");
			hs_JPsi->Draw("e hist same");
			Fit_func_MC.Draw_Functions();

			auto legend_acc = new TLegend(0.54, 0.87, 0.90, 0.60);
			legend_acc->AddEntry(h_only_JPsi, Form("nb JPsi %3.1f", nb_JPsi_integral), "f1");
			legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("J#psi fit (%3.1f #pm %3.1f)  ", Fit_func_MC.Get_Integral_Signal(), Fit_func_MC.Get_Integral_Error_Signal()), "l");
			legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("#Chi^{2} %3.1f, NdF %3.1f, #Chi^{2}/NdF %3.1f ", Fit_func_MC.chi2, Fit_func_MC.NDF, Fit_func_MC.chi2 / Fit_func_MC.NDF), "");
			legend_acc->SetFillStyle(0);
			legend_acc->SetLineWidth(0);
			legend_acc->Draw("same ");

			double nb_JPsi_MC = Fit_func_MC.Get_Integral_Signal(); // nb_JPsi_integral; //

			Acc_Num.push_back(nb_JPsi_MC); //(sample_Acc->Integral());
			// cout<<"acc num 1 "<<sample_Acc->Integral()<<endl;

			cancAcc->SaveAs(output_folder + name_pdf + ".pdf");

			*/
		}
	}

	void Process_MC_Diff()
	{

		for (int i = 0; i < labels_MC.size(); i++)
		{

			//////////////////////////////////////////////////
			// Set options for each label
			TString label = labels_MC[i][0];
			TString xAxis_label = labels_MC[i][1];
			TString min_histo_option = labels_MC[i][2];
			TString max_histo_option = labels_MC[i][3];
			TString nb_bins = labels_MC[i][4];
			TString string_cut = labels_MC[i][5];
			TString max_y = labels_MC[i][6];
			double variable_min = stof((string)labels_MC[i][7].Data());
			double variable_max = stof((string)labels_MC[i][8].Data());

			std::vector<double> integral_flux{};
			//////////////////////////////////////////////////

			cout << "\n";
			cout << "//////////////////////////////////////////////////"
				 << "\n";
			cout << "Doing " << label << " " << i << " plot"
				 << "\n";
			cout << "//////////////////////////////////////////////////"
				 << "\n";
			cout << "\n";

			TCut cut = string_cut.Data(); //""; //"1";////

			TString cut_string = cut.GetTitle();

			float min_histo = stof((string)min_histo_option.Data());
			float max_histo = stof((string)max_histo_option.Data());

			int nbBins = stoi((string)nb_bins.Data());

			if (debug)
				cout << "Bining and range " << nbBins << " " << min_histo << " " << max_histo << "\n";

			auto legend = new TLegend(0.54, 0.87, 0.90, 0.60);

			THStack *hs_MC = new THStack("hs_MC", "");
			THStack *hs_MC_no_rad = new THStack("hs_MC_no_rad", "");

			for (int j = 0; j < Analysis_Sample.samples.size(); j++)
			{

				double lumi_sample = stof((string)Analysis_Sample.samples[j][4].Data());
				double xsec = stof((string)Analysis_Sample.samples[j][1].Data());
				int nbEvents_sample = Analysis_Sample.ngen[j];

				if (Analysis_Sample.samples[j][3] == "J#psi")
				{

					if (debug)
						cout << "Nb of events in " << Analysis_Sample.samples[j][0] << " : " << nbEvents_sample << "\n";

					TString hist_name1 = Form("sample_hist1_%s_%i_%i", Analysis_Sample.samples[j][3].Data(), j, i);
					TH1D *sample_hist1 = new TH1D(hist_name1, "", nbBins, min_histo, max_histo);

					TString hist_Acc = Form("sample_Acc_%s_%i_%i", Analysis_Sample.samples[j][3].Data(), j, i);
					TH1D *sample_Acc = new TH1D(hist_Acc, "", nbBins, min_histo, max_histo);

					TCut weight1 = Form("(%s+%s)*%f", "virtual_flux_Gen", "real_flux_Gen", lumi_sample);
					TCut weight_acc = Form("%s*%f*%f*%f/(%i)", "weight", xsec, Analysis_Sample.lumi_factor, lumi_sample, nbEvents_sample);

					// MC_tree[j]->Draw(label + ">>" + hist_name, cut * weight);
					Analysis_Sample.MC_tree[j]->Draw(label + ">>" + hist_name1, cut * weight1 * weight_acc);
					Analysis_Sample.MC_tree[j]->Draw(label + ">>" + hist_Acc, cut * weight_acc);

					if (debug)
					{
						// cout << cut * weight << "\n";
						cout << "Debug flux"
							 << "\n";
						cout << weight1 << endl;
						// cout << "Number of entries " << sample_hist->GetEntries() << "\n";
						cout << "Number of entries " << sample_hist1->GetEntries() << "\n";
						cout << "Integral flux " << sample_hist1->Integral() << "\n";
						TCanvas *candebug = new TCanvas("", "candebug", 1500, 1000);
						sample_hist1->Draw();
						candebug->SaveAs(output_folder + name_pdf + ".pdf");
					}

					sample_Acc->SetLineWidth(0);
					sample_Acc->SetLineColor(kBlack);
					sample_Acc->SetMarkerSize(0);
					sample_Acc->SetFillColor(std::stoi(Analysis_Sample.samples[j][2].Data()));
					sample_Acc->SetStats(kFALSE);

					hs_MC->SetTitle(";" + label + ";Events");
					hs_MC->Add(sample_Acc);

					// double integral_flux_sample = (sample_hist1->GetEntries() > 100) ? (variable_max - variable_min) * ((sample_hist1->Integral()) / sample_hist1->GetEntries()) : 0.0;
					double integral_flux_sample = (sample_hist1->GetEntries() > 100) ? (Eg_max - Eg_min) * ((sample_hist1->Integral()) / sample_Acc->Integral()) : 0.0;
					// double integral_flux_sample = (sample_hist1->GetEntries() > 100) ? (Eg_max - Eg_min) * ((sample_hist1->Integral()) / sample_hist1->GetEntries()) : 0.0;
					integral_flux.push_back(integral_flux_sample);

					if (debug)
						cout << Analysis_Sample.samples[j][0] << "\n";

					legend->AddEntry(sample_Acc, Form("%s (%3.6f), #sigma=%3.1f pb", Analysis_Sample.samples[j][3].Data(), sample_Acc->Integral(), stof(Analysis_Sample.samples[j][1].Data())), "f1");
				}

				if (Analysis_Sample.samples[j][3] == "Jpsi_Rad_Corr")
				{
					TString hist_no_rad = Form("sample_no_rad_%s_%i_%i", Analysis_Sample.samples[j][3].Data(), j, i);
					TH1D *sample_no_rad = new TH1D(hist_no_rad, "", nbBins, min_histo, max_histo);

					TCut weight_no_rad = Form("%s*%f*%f*%f/(%i)", "weight", xsec, Analysis_Sample.lumi_factor, lumi_sample, nbEvents_sample);

					Analysis_Sample.MC_tree[j]->Draw(label + ">>" + hist_no_rad, cut * weight_no_rad);

					hs_MC_no_rad->SetTitle(";" + label + ";Events");
					hs_MC_no_rad->Add(sample_no_rad);
				}
			}

			/////////////////////////
			// Calculation of CS
			/////////////////////////(TH1D *)(hs->GetStack()->Last())
			// TH1D *flux_stack_hist = (TH1D *)(hs_MC_flux->GetStack()->Last());
			// double avg_flux_wrong = (Eg_max - Eg_min) * ((flux_stack_hist->Integral()) / flux_stack_hist->GetEntries()); //////////////////////This is right formula but need sample normalization
			// cout << "Average Flux old formula " << avg_flux_wrong << "\n";												 //// maybe use th1f of the flux and use getMean
			double avg_flux = std::accumulate(integral_flux.begin(), integral_flux.end(), 0.0);
			for (const auto &element : integral_flux)
			{
				std::cout << element << " ";
			}
			cout << "Average Flux " << avg_flux << "\n"; //// maybe use th1f of the flux and use getMean
			TH1D *MC_stack_hist = (TH1D *)(hs_MC->GetStack()->Last());
			double Acc = Acc_Num[i] / (Analysis_Sample.reduction_factor_Gen * MC_stack_hist->Integral());
			//double Acc = Acc_Num[i] / (MC_stack_hist->Integral());

			TH1D *No_rad_MC_stack_hist = (TH1D *)(hs_MC_no_rad->GetStack()->Last());
			double Rad_corr = (MC_stack_hist->Integral()) / (No_rad_MC_stack_hist->Integral());

			Acc_hist->SetBinContent(i + 1, Acc);
			Flux_hist->SetBinContent(i + 1, avg_flux * Analysis_Sample.lumi_factor);
			Nb_JPsi_hist->SetBinContent(i + 1, nb_JPsi_Data[i]);
			Rad_corr_hist->SetBinContent(i + 1, Rad_corr);

			cout << "Acc " << Acc << "\n";
			cout << "Nb JPsi " << nb_JPsi_Data[i] << "\n";
			cout << "Br " << Branching_ratio << "\n";
			cout << "w_c " << w_c[i] << "\n";
			cout << "Size of the bin Delta_variable " << (variable_max - variable_min) << "\n";
			cout << "Size of the bin E " << Eg_min << " " << Eg_max << " " << (Eg_max - Eg_min) << "\n";
			double bin_volume_corr = bin_volume_correction(variable_min, variable_max, Eg_min, Eg_max, debug, output_folder + name_pdf);
			cout << "Bin volume correction " << bin_volume_corr << endl;

			double CS_usual = 0.001 * nb_JPsi_Data[i] / ((variable_max - variable_min) * avg_flux * Analysis_Sample.lumi_factor * Branching_ratio * (Acc)*w_c[i] * bin_volume_corr);
			cout << "CS using normal formula " << CS_usual << "\n";						  // remove lumi factor ?
			double error_CS_usual = CS_usual * (error_nb_JPsi_Data[i] / nb_JPsi_Data[i]); // 0.001 * error_nb_JPsi_Data[i] / ((variable_max - variable_min) * avg_flux * Analysis_Sample.lumi_factor * Branching_ratio * (Acc)*w_c[i] * bin_volume_correction); //
			cout << "error on the jpsi number " << error_nb_JPsi_Data[i] << "\n";
			cout << "error CS using normal formula " << error_CS_usual << "\n";

			cout << "Nb JPsi (comptage)" << nb_JPsi_Data_C[i] << "\n";
			double CS_comptage = 0.001 * nb_JPsi_Data_C[i] / (avg_flux * Analysis_Sample.lumi_factor * Branching_ratio * (Acc)*w_c[i]);				//
			double error_CS_comptage = 0.001 * sqrt(nb_JPsi_Data_C[i]) / (avg_flux * Analysis_Sample.lumi_factor * Branching_ratio * (Acc)*w_c[i]); //
			cout << "CS using comptage formula " << CS_comptage << "\n";

			cout << "average variable " << average_variable[i] << "\n";

			Latex_Table.Add_value(i, labels_MC.size(), CS_usual, error_CS_usual, average_variable[i], sigma_variable[i], average_variable_2[i], sigma_variable_2[i]);

			JPsi_CS_Graph.SetPoint(i, average_variable[i], CS_usual);
			// JPsi_CS_Graph.SetPointError(i, average_variable[i] - variable_min, variable_max - average_variable[i], error_CS_usual, error_CS_usual);
			// JPsi_CS_Graph.SetPointError(i, sigma_variable[i], sigma_variable[i], error_CS_usual, error_CS_usual);
			JPsi_CS_Graph.SetPointError(i, 0.0, 0.0, error_CS_usual, error_CS_usual);

			JPsi_CS_Graph_C.SetPoint(i, average_variable[i], CS_usual);
			JPsi_CS_Graph_C.SetPointError(i, average_variable[i] - variable_min, variable_max - average_variable[i], error_CS_usual, error_CS_usual);
			/////////////////////////
			/////////////////////////

			TCanvas *cancG0 = new TCanvas("", "can0", 1500, 1000);

			float max_display = (hs_MC->GetMaximum()) * 1.5;
			hs_MC->SetMaximum(max_display);

			hs_MC->Draw("hist"); // nostack

			TText *t = new TText();
			t->SetTextAlign(22);
			t->SetTextColorAlpha(kGray, 0.50);
			t->SetTextFont(40);
			t->SetTextSize(0.25);
			t->SetTextAngle(25);

			legend->SetFillStyle(0);
			legend->SetLineWidth(0);
			legend->Draw("same ");

			// Labels
			double x_top_label = 0.90;

			TPaveText *CLAS12_Internal = new TPaveText(0.10, x_top_label, 0.288191, x_top_label + 0.1, "NDC");
			CLAS12_Internal->SetFillStyle(4050);
			CLAS12_Internal->SetLineColor(0);
			CLAS12_Internal->SetTextFont(42);
			CLAS12_Internal->SetTextSize(0.0599401);
			CLAS12_Internal->SetBorderSize(0);
			CLAS12_Internal->AddText(Form("CLAS12 Preliminary"));
			CLAS12_Internal->Draw();

			t->DrawTextNDC(.5, .53, "Preliminary");

			cancG0->SaveAs(output_folder + name_pdf + ".pdf");
		}
	}

	void Finalize_Diff_CS_Calculation()
	{
		TGraphAsymmErrors *JPsi_CS_GlueX;
		TGraphAsymmErrors *JPsi_CS_GlueX_Binned;

		if (bin_id == "bin 2")
		{
			JPsi_CS_GlueX = new TGraphAsymmErrors(9, t_GluX_Bin2, xs_t_GluX_Bin2, e_t_down_GluX_Bin2_zeros, e_t_up_GluX_Bin2_zeros, e_xs_t_GluX_Bin2, e_xs_t_GluX_Bin2);
			JPsi_CS_GlueX_Binned = new TGraphAsymmErrors(9, t_GluX_Bin2, xs_t_GluX_Bin2, e_t_down_GluX_Bin2, e_t_up_GluX_Bin2, e_xs_t_GluX_Bin2, e_xs_t_GluX_Bin2);
		}
		else
		{
			JPsi_CS_GlueX = new TGraphAsymmErrors(7, t_GluX_Bin1, xs_t_GluX_Bin1, e_t_down_GluX_Bin1_zeros, e_t_up_GluX_Bin1_zeros, e_xs_t_GluX_Bin1, e_xs_t_GluX_Bin1);
			JPsi_CS_GlueX_Binned = new TGraphAsymmErrors(7, t_GluX_Bin1, xs_t_GluX_Bin1, e_t_down_GluX_Bin1, e_t_up_GluX_Bin1, e_xs_t_GluX_Bin1, e_xs_t_GluX_Bin1);
		}

		JPsi_CS_GlueX->SetMarkerColor(kRed);
		JPsi_CS_GlueX->SetLineColor(kRed);
		JPsi_CS_GlueX->SetMarkerStyle(22);
		JPsi_CS_GlueX->SetMarkerSize(3.);
		JPsi_CS_GlueX_Binned->SetMarkerColor(kRed);
		JPsi_CS_GlueX_Binned->SetLineColor(kRed);
		JPsi_CS_GlueX_Binned->SetMarkerStyle(22);
		JPsi_CS_GlueX_Binned->SetMarkerSize(3.);

		/*TGraphAsymmErrors *JPsi_CS_GlueX = new TGraphAsymmErrors(7, t_GluX_Bin1, xs_t_GluX_Bin1, e_t_down_GluX_Bin1, e_t_up_GluX_Bin1, e_xs_t_GluX_Bin1, e_xs_t_GluX_Bin1);
		JPsi_CS_GlueX->SetMarkerColor(kRed);
		JPsi_CS_GlueX->SetLineColor(kRed);
		JPsi_CS_GlueX->SetMarkerStyle(22);
		JPsi_CS_GlueX->SetMarkerSize(3.);*/

		/*TGraphAsymmErrors *JPsi_CS_GlueX = new TGraphAsymmErrors(7, t_GluX_Bin1, xs_t_GluX_Bin1, e_t_down_GluX_Bin1_zeros, e_t_up_GluX_Bin1_zeros, e_xs_t_GluX_Bin1, e_xs_t_GluX_Bin1);
		JPsi_CS_GlueX->SetMarkerColor(kRed);
		JPsi_CS_GlueX->SetLineColor(kRed);
		JPsi_CS_GlueX->SetMarkerStyle(22);
		JPsi_CS_GlueX->SetMarkerSize(3.);*/

		TF1 *f_JPsi_sigm_int_2g = new TF1("f_JPsi_sigm_int_2g", JPsi_sigm_int_2g, 8.25, 10.7, 3);
		double SLAC_Fit_scale = 7.79117e-23;
		double tSlope = 1.13;
		f_JPsi_sigm_int_2g->SetParameter(1, SLAC_Fit_scale);
		f_JPsi_sigm_int_2g->SetParameter(2, tSlope);

		/*TFile *CS_file = new TFile("../Spring2019/outputTCS_JPsi_Spring2019_Pass2.root");
		TTree *CS_tree = (TTree *)CS_file->Get("tree_Gen");
		TH1F *CS_hist_1 = new TH1F("CS_hist_1", "CS_hist_1", 50, 8.6, 10.2);
		CS_tree->Draw("Epho_Gen>>CS_hist_1", "weight*(1.0)/(((10.2-8.6)/50.)*1124999*(virtual_flux_Gen+real_flux_Gen))", "hist");*/

		cout << "here" << endl;

		cout << "here" << endl;
		///////////////////////////////////////////////////////////////////////////////

		TCanvas *cancG0 = new TCanvas("", "can0", 1500, 1000);
		JPsi_CS_Graph.SetTitle(";-t [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");
		JPsi_CS_Graph.GetXaxis()->SetRangeUser(0.0, 8.0);
		JPsi_CS_Graph.SetMarkerColor(kBlue);
		JPsi_CS_Graph.SetLineColor(kBlue);
		JPsi_CS_Graph.SetMarkerStyle(21);
		JPsi_CS_Graph.SetMarkerSize(2.);
		JPsi_CS_Graph_C.SetMarkerColor(kOrange);
		JPsi_CS_Graph_C.SetLineColor(kOrange);
		JPsi_CS_Graph_C.SetMarkerStyle(21);
		JPsi_CS_Graph_C.SetMarkerSize(2.);
		// JPsi_CS_Graph.SetMinimum(0.0005);
		JPsi_CS_Graph.SetMaximum(3.0);
		// JPsi_CS_Graph.SetMaximum(1.);
		// JPsi_CS_Graph.SetMinimum(0.);
		// JPsi_CS_Graph.SetMaximum(1.8);
		// JPsi_CS_GlueX->Draw("AP");
		// JPsi_CS_HallC->Draw("P");
		// JPsi_CS_GlueX_old->Draw("P");
		gPad->SetLogy();
		JPsi_CS_Graph.Draw("AP");
		auto legend_CS = new TLegend(0.40, 0.87, 0.90, 0.60);
		legend_CS->AddEntry(&JPsi_CS_Graph, "CLAS12", "lp");
		// legend_CS->AddEntry(JPsi_CS_HallC, "Hall C (unpublished)", "lp");
		// legend_CS->AddEntry(JPsi_CS_GlueX_old, "GlueX (2019)", "lp");
		legend_CS->SetFillStyle(0);
		legend_CS->SetLineWidth(0);
		legend_CS->Draw("same");

		TText *t = new TText();
		t->SetTextAlign(22);
		t->SetTextColorAlpha(kGray, 0.70);
		t->SetTextFont(40);
		t->SetTextSize(0.25);
		t->SetTextAngle(25);
		t->DrawTextNDC(.5, .53, "Preliminary");

		cancG0->SaveAs(output_folder + name_pdf + ".pdf");

		/// Log CS
		TCanvas *cancG01 = new TCanvas("", "can01", 1500, 1000);
		JPsi_CS_GlueX->SetTitle(";-t [GeV^{2}]; d#sigma/dt [nb/GeV^{2}]");
		gPad->SetLogy();

		TF1 *func_fit_t = new TF1("func_fit_t", Fit_diff_CS_t, 0, 7, 2);
		func_fit_t->SetParameters(1.0, 1.5);
		func_fit_t->SetLineColor(kCyan);

		double max_fit_glueX = 3.5;

		if (bin_id == "bin 2")
			max_fit_glueX = 6.5;

		TF1 *func_fit_t_2 = new TF1("func_fit_t_2", Fit_diff_CS_t, 0, max_fit_glueX, 2);
		func_fit_t_2->SetParameters(1.0, 1.5);
		func_fit_t_2->SetLineColor(kOrange);

		JPsi_CS_GlueX->SetMaximum(3.0);
		JPsi_CS_GlueX->SetMinimum(0.001);
		JPsi_CS_GlueX->Draw("AP");
		JPsi_CS_GlueX->Fit(func_fit_t_2, "R");

		JPsi_CS_Graph.Draw("P");
		JPsi_CS_Graph.Fit(func_fit_t);

		double sigma_0_CLAS12 = func_fit_t->GetParameter(0);
		double m_s_CLAS12 = func_fit_t->GetParameter(1);
		double e_sigma_0_CLAS12 = func_fit_t->GetParError(0);
		double e_m_s_CLAS12 = func_fit_t->GetParError(1);
		double r_m_CLAS12 = convert_Gev_to_fm * sqrt(12.) / m_s_CLAS12;
		double e_r_m_CLAS12 = r_m_CLAS12 * e_m_s_CLAS12 / m_s_CLAS12;

		double sigma_0_GlueX = func_fit_t_2->GetParameter(0);
		double m_s_GlueX = func_fit_t_2->GetParameter(1);
		double e_sigma_0_GlueX = func_fit_t_2->GetParError(0);
		double e_m_s_GlueX = func_fit_t_2->GetParError(1);
		double r_m_GlueX = convert_Gev_to_fm * sqrt(12.) / m_s_GlueX;
		double e_r_m_GlueX = r_m_GlueX * e_m_s_GlueX / m_s_GlueX;

		legend_CS->AddEntry(JPsi_CS_GlueX, "GlueX (2023)", "lp");
		legend_CS->AddEntry(func_fit_t, Form("Fit CLAS12 (#sigma_{0} %3.2f #pm %3.2f, m_{s} %3.2f #pm %3.2f GeV)  ", sigma_0_CLAS12, e_sigma_0_CLAS12, m_s_CLAS12, e_m_s_CLAS12), "l");
		legend_CS->AddEntry(func_fit_t, Form("r_{m} %3.2f #pm %3.2f fm", r_m_CLAS12, e_r_m_CLAS12), "");
		legend_CS->AddEntry(func_fit_t_2, Form("Fit GlueX (#sigma_{0} %3.2f #pm %3.2f, m_{s} %3.2f #pm %3.2f GeV)  ", sigma_0_GlueX, e_sigma_0_GlueX, m_s_GlueX, e_m_s_GlueX), "l");
		legend_CS->AddEntry(func_fit_t_2, Form("r_{m} %3.2f #pm %3.2f fm", r_m_GlueX, e_r_m_GlueX), "");
		legend_CS->Draw("same");

		t->DrawTextNDC(.5, .53, "Preliminary");

		cancG01->SaveAs(output_folder + name_pdf + ".pdf");

		/// Integrated CS/////////////
		TCanvas *cancG1 = new TCanvas("", "can1", 1500, 1000);

		TGraphErrors *JPsi_CS_GlueX_Eg = new TGraphErrors(18, EgGluex, xsGluex, xeGluex, sxsGluex);
		JPsi_CS_GlueX_Eg->SetMarkerColor(kRed);
		JPsi_CS_GlueX_Eg->SetLineColor(kRed);
		JPsi_CS_GlueX_Eg->SetMarkerStyle(22);
		JPsi_CS_GlueX_Eg->SetMarkerSize(3.);
		JPsi_CS_GlueX_Eg->SetTitle(";Eg [GeV]; #sigma [nb]");

		double integrated_CS_CLAS12 = 0.0;
		double error_integrated_CS_CLAS12 = 0.0;

		double integrated_CS_GLUEX = 0.0;
		double error_integrated_CS_GLUEX = 0.0;

		auto JPsi_CS_Integrated = new TGraphErrors(1);
		auto JPsi_CS_Integrated_GLUEX = new TGraphErrors(1);

		JPsi_CS_Integrated->SetMarkerColor(kCyan);
		JPsi_CS_Integrated->SetLineColor(kCyan);
		JPsi_CS_Integrated->SetMarkerStyle(21);
		JPsi_CS_Integrated->SetMarkerSize(2.);

		JPsi_CS_Integrated_GLUEX->SetMarkerColor(kOrange);
		JPsi_CS_Integrated_GLUEX->SetLineColor(kOrange);
		JPsi_CS_Integrated_GLUEX->SetMarkerStyle(21);
		JPsi_CS_Integrated_GLUEX->SetMarkerSize(2.);

		for (int i = 0; i < JPsi_CS_Graph.GetN(); i++)
		{
			integrated_CS_CLAS12 += JPsi_CS_Graph_C.GetPointY(i) * (JPsi_CS_Graph_C.GetErrorXhigh(i) + JPsi_CS_Graph_C.GetErrorXlow(i));
			cout << "bin " << i << " width " << (JPsi_CS_Graph_C.GetErrorXhigh(i) + JPsi_CS_Graph_C.GetErrorXlow(i)) << endl;
			double error_contribution = JPsi_CS_Graph_C.GetErrorYhigh(i) * (JPsi_CS_Graph_C.GetErrorXhigh(i) + JPsi_CS_Graph_C.GetErrorXlow(i));
			error_integrated_CS_CLAS12 += error_contribution * error_contribution;
		}
		JPsi_CS_Integrated->SetPoint(0, Eg_min + (Eg_max - Eg_min) / 2.0, integrated_CS_CLAS12);
		JPsi_CS_Integrated->SetPointError(0, (Eg_max - Eg_min) / 2.0, sqrt(error_integrated_CS_CLAS12));
		cout << "integrated CS " << integrated_CS_CLAS12 << endl;

		for (int i = 0; i < JPsi_CS_GlueX_Binned->GetN(); i++)
		{
			integrated_CS_GLUEX += JPsi_CS_GlueX_Binned->GetPointY(i) * (JPsi_CS_GlueX_Binned->GetErrorXhigh(i) + JPsi_CS_GlueX_Binned->GetErrorXlow(i));
			cout << "bin " << i << " width " << (JPsi_CS_GlueX_Binned->GetErrorXhigh(i) + JPsi_CS_GlueX_Binned->GetErrorXlow(i)) << endl;
			double error_contribution = JPsi_CS_GlueX_Binned->GetErrorYhigh(i) * (JPsi_CS_GlueX_Binned->GetErrorXhigh(i) + JPsi_CS_GlueX_Binned->GetErrorXlow(i));
			error_integrated_CS_GLUEX += error_contribution * error_contribution;
		}
		JPsi_CS_Integrated_GLUEX->SetPoint(0, Eg_min + (Eg_max - Eg_min) / 2.0, integrated_CS_GLUEX);
		JPsi_CS_Integrated_GLUEX->SetPointError(0, (Eg_max - Eg_min) / 2.0, error_integrated_CS_GLUEX);
		cout << "integrated CS " << integrated_CS_GLUEX << endl;

		TFile *CS_CLAS12 = new TFile("../CS_Extraction_combine_good_error_Gaussian_with_int_SLER_CS_graph.root");
		TGraphAsymmErrors *CS_CLAS12_Graph = (TGraphAsymmErrors *)CS_CLAS12->Get("Graph");

		JPsi_CS_GlueX_Eg->Draw("AP");
		JPsi_CS_Integrated->Draw("P");
		JPsi_CS_Integrated_GLUEX->Draw("P");
		CS_CLAS12_Graph->Draw("P");
		// f_JPsi_sigm_int_2g->Draw("same");
		// CS_graph_theo->Draw("same");
		auto legend_CS_1 = new TLegend(0.12, 0.7, 0.48, 0.85);
		legend_CS_1->AddEntry(JPsi_CS_Integrated, "Integrated CS CLAS12", "lp");
		legend_CS_1->AddEntry(JPsi_CS_Integrated_GLUEX, "Integrated CS GLUEX", "lp");
		// legend_CS->AddEntry(JPsi_CS_HallC, "Hall C (unpublished)", "lp");
		// legend_CS_1->AddEntry(JPsi_CS_GlueX, "GlueX (2023)", "lp");
		// legend_CS_1->AddEntry(JPsi_CS_GlueX_old, "GlueX (2019)", "lp");
		legend_CS_1->SetFillStyle(0);
		legend_CS_1->SetLineWidth(0);
		legend_CS_1->Draw("same");

		t->DrawTextNDC(.5, .53, "Preliminary");

		cancG1->SaveAs(output_folder + name_pdf + ".pdf");

		gStyle->SetPaintTextFormat("4.3f");
		TCanvas *cancG2 = new TCanvas("", "can2", 1500, 1000);
		Acc_hist->SetStats(kFALSE);
		Acc_hist->GetYaxis()->SetTitleOffset(1.3);
		Acc_hist->SetTitle(";Bin;Acc");
		Acc_hist->Draw("hist text0");
		cancG2->SaveAs(output_folder + name_pdf + ".pdf");
		TCanvas *cancG3 = new TCanvas("", "can3", 1500, 1000);
		Flux_hist->SetStats(kFALSE);
		Flux_hist->GetYaxis()->SetTitleOffset(1.3);
		Flux_hist->SetTitle(";Bin;Flux");
		Flux_hist->Draw("hist text0");
		cancG3->SaveAs(output_folder + name_pdf + ".pdf");
		TCanvas *cancG5 = new TCanvas("", "can5", 1500, 1000);
		Nb_JPsi_hist->SetStats(kFALSE);
		Nb_JPsi_hist->GetYaxis()->SetTitleOffset(1.3);
		Nb_JPsi_hist->SetTitle(";Bin;Nb JPsi");
		Nb_JPsi_hist->Draw("hist text0");
		cancG5->SaveAs(output_folder + name_pdf + ".pdf");
		TCanvas *cancG4 = new TCanvas("", "can4", 1500, 1000);
		Rad_corr_hist->SetStats(kFALSE);
		Rad_corr_hist->GetYaxis()->SetTitleOffset(1.3);
		Rad_corr_hist->SetTitle(";Bin;Rad. Corr.");
		Rad_corr_hist->Draw("hist text0");
		cancG4->SaveAs(output_folder + name_pdf + ".pdf)");
	}
	//////////////////////////////
};

#endif
