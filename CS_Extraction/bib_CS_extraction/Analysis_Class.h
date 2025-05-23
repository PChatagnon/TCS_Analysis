#ifndef Analysis_Class
#define Analysis_Class

#include "Utils.h"

class Analysis
{
public: // Keep everything public for convenience
	////// Store configuration //////
	TString RunGroup = "RGA";
	TString variable = "Epho";
	TString bin_id = "bin 1";
	TString fit_procedure = "Default";
	bool debug = false;
	bool ratio_pad = false;
	bool pull_pad = true;
	/////////////////////////////////

	////// Store fit and plotting configuration //////
	double min_fit = 2.7;
	double max_fit = 3.3;
	/////////////////////////////////

	////// Store normalization configuration //////
	double Mass_norm_low = 2.6;
	double Mass_norm_high = 2.9;
	/////////////////////////////////

	////// Store all quantities for the analysis //////
	std::vector<vector<TString>> labels{};
	std::vector<vector<TString>> labels_MC{};

	std::vector<double> nb_JPsi_Data{};
	std::vector<double> error_nb_JPsi_Data{};

	std::vector<double> nb_JPsi_Data_C{};
	std::vector<double> error_nb_JPsi_Data_C{};

	std::vector<double> w_c{};
	std::vector<double> Acc_Num{};
	std::vector<double> average_variable{};
	std::vector<double> sigma_variable{};
	/////////////////////////////////

	////// Store the path to data and output //////
	TString output_folder = "/Users/pc281174/Desktop/JPsi_analysis/CS_Extraction/";
	TString name_pdf = "CS_Extraction_Gauss_exp_1";
	/////////////////////////////////

	////// Store kinematic and exclusivity cuts //////
	TCut kinematic_cut = "pass_EC_cut &&  Proton.Theta()*180./3.141592<35. && M>2.6 && (Electron.P() > 1.7) && (Positron.P() > 1.7) && positron_SF>0.15 && electron_SF>0.15 && ( Positron.P()<4.0 || (Positron.P()>4.0 && positron_score>0.00)) && ( Electron.P()<4.0 || (Electron.P()>4.0 && electron_score>0.00))  && positron_HTCC_ECAL_match==1. && electron_HTCC_ECAL_match==1.";
	TCut kinematic_cut_BG = "Proton.Theta()*180./3.141592<35. &&  M>2.6 && (Electron.P() > 1.7) && (Positron.P() > 1.7)";
	TCut exclusivity_cut = "abs(MMassBeam)<0.4  && abs(Q2)<0.5";
	TCut data_cut = "abs(positron_HTCCt-electron_HTCCt)<4";
	/////////////////////////////////

	////// Store histogramms and Graphs //////
	TH1D *Acc_hist;
	TH1D *Flux_hist;
	TH1D *Nb_JPsi_hist;
	TH1D *Rad_corr_hist;
	TGraphAsymmErrors JPsi_CS_Graph;
	/////////////////////////////////

	/////// Store the Latex table //////
	Latex_Table_writter Latex_Table;

	////// Store dataset //////
	Sample Analysis_Sample;

	//// Default creator
	Analysis() {}

	/////// Set configurations //////
	void Set_Sample_to_RGA()
	{
		Analysis_Sample.Setup_RGA();
		RunGroup = "RGA";
	}

	void Set_Sample_to_RGB()
	{
		Analysis_Sample.Setup_RGB();
		RunGroup = "RGB";
	}

	void Set_name_pdf(TString name_pdf_in)
	{
		name_pdf = name_pdf_in;
	}
	/////// End set configurations //////

	/////// Set parameters //////
	std::map<std::string, std::string> parseInputFile(const std::string &filename)
	{
		std::map<std::string, std::string> parameters;
		std::ifstream infile(filename);
		std::string line;

		while (std::getline(infile, line))
		{
			std::istringstream iss(line);
			std::string key, value;
			if (!(iss >> key))
			{
				// Skip empty lines
				continue;
			}

			// Check if the next character is a double quote
			char c;
			iss >> std::ws >> c;
			if (c != '"')
			{
				// If the value is not quoted, skip to the next line
				continue;
			}

			// Retrieve the value within double quotes
			std::getline(iss >> std::ws, value, '"');
			parameters[key] = value;
		}

		return parameters;
	}

	bool isParameter(const std::map<std::string, std::string> &parameters, const std::string &key)
	{
		auto it = parameters.find(key);
		return (it != parameters.end());
	}

	std::string getParameterValue(const std::map<std::string, std::string> &parameters, const std::string &key)
	{
		auto it = parameters.find(key);
		if (it != parameters.end())
		{
			return it->second;
		}
		else
		{
			return ""; // Return empty string if key is not found
		}
	}

	void Set_Parameters(string parameters_file)
	{
		std::map<std::string, std::string> parameters = parseInputFile(parameters_file);

		if (isParameter(parameters, "output_folder"))
			output_folder = TString(getParameterValue(parameters, "output_folder"));

		if (isParameter(parameters, "name_pdf"))
			name_pdf = TString(getParameterValue(parameters, "name_pdf"));

		if (isParameter(parameters, "fit_procedure"))
			fit_procedure = TString(getParameterValue(parameters, "fit_procedure"));

		if (isParameter(parameters, "bin_id"))
			bin_id = TString(getParameterValue(parameters, "bin_id"));

		if (isParameter(parameters, "exclusivity_cut"))
			exclusivity_cut = Form("%s", getParameterValue(parameters, "exclusivity_cut").c_str());

		if (isParameter(parameters, "kinematic_cut"))
			kinematic_cut = Form("%s", getParameterValue(parameters, "kinematic_cut").c_str());

		if (isParameter(parameters, "kinematic_cut_BG"))
			kinematic_cut_BG = Form("%s", getParameterValue(parameters, "kinematic_cut_BG").c_str());

		if (isParameter(parameters, "data_cut"))
			data_cut = Form("%s", getParameterValue(parameters, "data_cut").c_str());

		if (isParameter(parameters, "min_fit"))
			min_fit = stof(getParameterValue(parameters, "min_fit"));

		if (isParameter(parameters, "max_fit"))
			max_fit = stof(getParameterValue(parameters, "max_fit"));

		if (isParameter(parameters, "debug"))
			debug = (getParameterValue(parameters, "debug") == "true");

		if (isParameter(parameters, "pull_pad"))
			pull_pad = (getParameterValue(parameters, "pull_pad") == "true");
	}

	//////  End Set parameters  //////

	//////  Setup binnings  //////
	void
	Set_Binning_RGA()
	{

		// Store the binning and plotting configuration
		TString min_hist_MC = "1.8";
		TString min_hist = "2.6";
		TString max_hist = "3.5";
		TString bin_hist = "30.";

		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>8.2 && Epho<8.65 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>8.65 && Epho<8.9 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>8.9 && Epho<9.05 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.05 && Epho<9.2 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.2 && Epho<9.46 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.46 && Epho<9.7 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.7 && Epho<10. ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>10. && Epho<10.2 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>10.2 && Epho<10.4 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>10.4 && Epho<10.6 ", "", "M2"});

		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>8.2 && Epho_Gen<8.65 ", "50", "8.2", "8.65"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>8.65 && Epho_Gen<8.9 ", "50", "8.65", "8.9"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>8.9 && Epho_Gen<9.05 ", "50", "8.9", "9.05"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.05 && Epho_Gen<9.2 ", "50", "9.05", "9.2"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.2 && Epho_Gen<9.46", "50", "9.2", "9.46"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.46 && Epho_Gen<9.7 ", "50", "9.46", "9.7"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.7 && Epho_Gen<10.0 ", "50", "9.7", "10."});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>10.0 && Epho_Gen<10.2 ", "50", "10.", "10.2"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>10.2 && Epho_Gen<10.4 ", "50", "10.2", "10.4"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>10.4 && Epho_Gen<10.6 ", "50", "10.4", "10.6"});
	}

	void Set_Binning_RGB()
	{

		// Store the binning and plotting configuration
		TString min_hist_MC = "1.8";
		TString min_hist = "2.6";
		TString max_hist = "3.5";
		TString bin_hist = "30.";

		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>8.2 && Epho<9.0 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.0 && Epho<9.2 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.2 && Epho<9.4 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.4 && Epho<9.6 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.6 && Epho<9.8 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>9.8 && Epho<10.0 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>10.0 && Epho<10.2 ", "", "M2"});
		labels.push_back({"M", "M_{ee}", min_hist, max_hist, bin_hist, "M>2.0 && Epho>10.2 && Epho<10.6 ", "", "M2"});

		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>8.2 && Epho_Gen<9.0 ", "50", "8.2", "9.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.0 && Epho_Gen<9.2 ", "50", "9.0", "9.2"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.2 && Epho_Gen<9.4 ", "50", "9.2", "9.4"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.4 && Epho_Gen<9.6 ", "50", "9.4", "9.6"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.6 && Epho_Gen<9.8", "50", "9.6", "9.8"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>9.8 && Epho_Gen<10.0 ", "50", "9.8", "10.0"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>10.0 && Epho_Gen<10.2 ", "50", "10.0", "10.2"});
		labels_MC.push_back({"M_Gen_2", "M_{ee}", min_hist_MC, max_hist, bin_hist, "weight<100 && M_Gen_2>2.0 && Epho_Gen>10.2 && Epho_Gen<10.6 ", "50", "10.2", "10.6"});
	}
	//////////////////////////////

	////// Method histogramms ///////
	void Setup_Histo()
	{
		JPsi_CS_Graph.Set(labels.size());

		Acc_hist = new TH1D("Acc_hist", "Acc_hist", labels.size(), 0.0, labels.size());
		Acc_hist->SetStats(kFALSE);
		Acc_hist->GetYaxis()->SetTitleOffset(1.3);
		Acc_hist->SetTitle(";Bin;Acc");

		Flux_hist = new TH1D("Flux_hist", "Flux_hist", labels.size(), 0.0, labels.size());
		Flux_hist->SetStats(kFALSE);
		Flux_hist->GetYaxis()->SetTitleOffset(1.3);
		Flux_hist->SetTitle(";Bin;Flux");

		Nb_JPsi_hist = new TH1D("Nb_JPsi_hist", "Nb_JPsi_hist", labels.size(), 0.0, labels.size());
		Nb_JPsi_hist->SetStats(kFALSE);
		Nb_JPsi_hist->GetYaxis()->SetTitleOffset(1.3);
		Nb_JPsi_hist->SetTitle(";Bin;Nb JPsi");

		Rad_corr_hist = new TH1D("Rad_corr_hist", "Rad_corr_hist", labels.size(), 0.0, labels.size());
		Rad_corr_hist->SetStats(kFALSE);
		Rad_corr_hist->GetYaxis()->SetTitleOffset(1.3);
		Rad_corr_hist->SetTitle(";Bin;Rad. Corr.");

		JPsi_CS_Graph.SetTitle(";E_{#gamma} [GeV]; #sigma [nb]");
		JPsi_CS_Graph.SetMarkerColor(kBlue);
		JPsi_CS_Graph.SetLineColor(kBlue);
		JPsi_CS_Graph.SetMarkerStyle(21);
		JPsi_CS_Graph.SetMarkerSize(2.);
		JPsi_CS_Graph.SetMinimum(0.);
		JPsi_CS_Graph.SetMaximum(2.0);
	}

	void Save_Histo_to_root()
	{
		JPsi_CS_Graph.SaveAs(output_folder + name_pdf + "_CS_graph.root");
		Acc_hist->SaveAs(output_folder + name_pdf + "_Acc.root");
		Flux_hist->SaveAs(output_folder + name_pdf + "_Flux.root");
		Nb_JPsi_hist->SaveAs(output_folder + name_pdf + "_Nb_JPsi.root");
		Rad_corr_hist->SaveAs(output_folder + name_pdf + "_Rad_corr.root");
	}
	//////////////////////////////

	////// Method for Latex Tables //////
	void Setup_Latex_Table()
	{
		Latex_Table = Latex_Table_writter("", output_folder.Data(), "Epho");
		Latex_Table.Set_output_name(((string)name_pdf.Data()) + "_Latex_Table.txt");
	}
	//////////////////////////////

	void redefineErrors(TH1D *hist)
	{
		if (!hist)
			return; // Check if the histogram is valid

		// Loop over all bins (from 1 to Nbins, 0 is underflow, N+1 is overflow)
		for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
		{
			double content = hist->GetBinContent(bin);			   // Get the bin content
			double error = (content > 0) ? std::sqrt(content) : 0; // Calculate new error as sqrt(content)
			hist->SetBinError(bin, error);						   // Set the new error for the bin
		}
	}

	/////// Processing function for the integrated CS //////
	void Process_REC()
	{

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
			cout << "Doing " << i << " " << label << " plot"
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

			TH1D *Data_hist = new TH1D(Form("Data_hist_%i", i), "Data_hist", nBins, min_histo_ini, max_histo_ini);
			TH1D *Average_variable = new TH1D(Form("Average_variable_%i", i), "Average_variable", 100, 0., 11.);

			TCut weight_data = Form("%s", "weight");

			Analysis_Sample.filtered_Data_tree->Draw(label + ">>" + Form("Data_hist_%i", i), weight_data * data_cut * cut * exclusivity_cut * kinematic_cut);
			Analysis_Sample.filtered_Data_tree->Draw(variable + ">>" + Form("Average_variable_%i", i), weight_data * data_cut * cut * exclusivity_cut * kinematic_cut);
			average_variable.push_back(Average_variable->GetMean());
			sigma_variable.push_back(Average_variable->GetRMS());

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

			TH1D *BG_hist = new TH1D(Form("BG_hist_%i", i), "BG_hist", nBins, min_histo, max_histo);
			TH1D *BG_only_hist = new TH1D(Form("BG_only_hist_%i", i), "BG_only_hist", nBins, min_histo, max_histo);

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
				TCut weight = Form("%s*%f*%f*%f/(%i)", "weight", xsec, Analysis_Sample.lumi_factor, lumi_sample, nbEvents_sample); // Form("%f",1.0);//

				Analysis_Sample.reduced_samples_tree[j]->Draw(label + ">>" + hist_name, cut * weight);

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
				sample_hist->SetFillColor(std::stoi(Analysis_Sample.samples[j][2].Data()));
				// sample_hist->SetFillColorAlpha(std::stoi(Analysis_Sample.samples[j][2].Data()), 0.50);
				sample_hist->SetStats(kFALSE);

				hs->SetTitle(";" + label + ";Events");
				hs->Add(sample_hist);

				hs_JPsi->SetTitle(";" + label + ";Events");
				if (Analysis_Sample.samples[j][3] == "J#psi")
					hs_JPsi->Add(sample_hist);

				if (debug)
					cout << Analysis_Sample.samples[j][0] << "\n";
				BG_hist->Add(sample_hist);

				if (Analysis_Sample.samples[j][3] == "BG")
					BG_only_hist->Add(sample_hist);
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
			TPad *pad1 = new TPad(Form("pad1%i", i), Form("pad1%i", i), 0, limit_lower_pad, 1, 1.0);

			if (ratio_pad || pull_pad)
				pad1->SetBottomMargin(0.003);

			float max_display = (hs->GetMaximum()) * 1.5;
			hs->SetMaximum(max_display);

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
				Data_hist->SetMaximum(Data_hist->GetMaximum() * 1.3);

			Data_hist->Draw("e");

			cout << "///////////////////////" << endl;
			cout << "FIT DATA" << endl;
			cout << "///////////////////////" << endl;
			Fit_Function Fit_func;
			Fit_func.Set_Data_hist(Data_hist);
			Fit_func.Set_Limits(min_fit, max_fit);
			// Fit_func.Single_Gaussian_fit("SLER", Form("func_%i", i));
			// Fit_func.Single_Gaussian_Int_fit("SLER", Form("func_%i", i));

			if (fit_procedure == "Default")
				Fit_func.Single_Gaussian_Int_fit("SLER", Form("func_%i", i));

			if (fit_procedure == "Crystall ball Pol 2 BG")
				Fit_func.Crystall_Ball_fit("SLER", Form("func_%i", i));

			if (fit_procedure == "Crystall ball exp BG")
				Fit_func.Crystall_Ball_fit_exp("SLER", Form("func_%i", i));

			if (fit_procedure == "Pol 2 BG")
				Fit_func.Single_Gaussian_Int_fit_Pol_BG_V2("SLER", Form("func_%i", i));

			if (fit_procedure == "Double Gaussian")
				Fit_func.Double_Gaussian_Fit("SLR", Form("func_%i", i));
			// Fit_func.Single_Gaussian_Fit_Flat_BG("SLR",Form("func_%i", i));
			double chi2 = Fit_func.chi2;
			double NDF = Fit_func.NDF;
			cout << "///////////////////////" << endl;

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

			double w_c_norm_wrong = (Data_hist->Integral(Data_hist->FindBin(Mass_norm_low), Data_hist->FindBin(Mass_norm_high)) / BG_hist->Integral(BG_hist->FindBin(Mass_norm_low), BG_hist->FindBin(Mass_norm_high)));
			cout << "W_C wrong " << w_c_norm_wrong << endl;

			double nb_data = (Data_hist->Integral(Data_hist->FindBin(Mass_norm_low), Data_hist->FindBin(Mass_norm_high)));
			double nb_bg_only = (BG_only_hist->Integral(BG_only_hist->FindBin(Mass_norm_low), BG_only_hist->FindBin(Mass_norm_high)));
			double nb_bg = (BG_hist->Integral(BG_hist->FindBin(Mass_norm_low), BG_hist->FindBin(Mass_norm_high)));
			double w_c_norm = (nb_data - nb_bg_only) / (nb_bg - nb_bg_only);
			cout << "W_C true ingredient nb_data " << nb_data << " nb_bg_only " << nb_bg_only << " nb_bg " << nb_bg << endl;
			cout << "W_C true " << w_c_norm << endl;
			w_c.push_back(w_c_from_BG_estimation); // 1.0); // 0.7); //.push_back(w_c_norm); // Data_hist->Integral(Data_hist->FindBin(2.6), Data_hist->FindBin(2.9)) / BG_hist->Integral(BG_hist->FindBin(2.6), BG_hist->FindBin(2.9)));

			legend->AddEntry(Fit_func.function_Signal, Form("J#psi fit (%3.1f #pm %3.1f)  ", nb_JPsi, error_nb_JPsi), "l");
			legend->AddEntry(Fit_func.function_Signal, Form("#Chi^{2}: %3.1f, NdF: %3.1f, #Chi^{2}/NdF: %3.1f ", chi2, NDF, chi2 / NDF), "");

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

			// lower plot will be in pad
			///////////////////////////////////////////////////////////////////////////////
			if (ratio_pad)
			{
				cancG0->cd(); // Go back to the main canvas before defining pad2
				TPad *pad2 = new TPad(Form("pad2%i", i), Form("pad2%i", i), 0, 0.05, 1, limit_lower_pad);
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
				TPad *pad2 = new TPad(Form("pad3%i", i), Form("pad3%i", i), 0, 0.05, 1, limit_lower_pad);
				pad2->SetTopMargin(0);
				pad2->SetGridy();
				pad2->SetBottomMargin(0.3);
				pad2->Draw();
				pad2->cd();

				TH1F *pulls_histo = new TH1F(Form("pullsHistogram%i", i), "", Data_hist->GetNbinsX(), Data_hist->GetXaxis()->GetXmin(), Data_hist->GetXaxis()->GetXmax());
				for (int i = 1; i <= Data_hist->GetNbinsX(); ++i)
				{
					double dataValue = Data_hist->GetBinContent(i);
					double dataError = Data_hist->GetBinError(i);
					double fitValue = Fit_func.function->Eval(Data_hist->GetBinCenter(i));
					double pull = (dataValue - fitValue) / dataError;
					if (dataError == 0.0)
						pull = 0.0;
					// cout << "pull " << pull << endl;
					pulls_histo->SetBinContent(i, pull);
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
				// pulls_histo->GetXaxis()->SetTitleOffset(4.);
				pulls_histo->GetXaxis()->SetLabelFont(43);
				pulls_histo->GetXaxis()->SetLabelSize(30);
			}
			///////////////////////////////////////////////////////////////////////////////

			if (i == 0)
				cancG0->SaveAs(output_folder + name_pdf + ".pdf(");
			else
				cancG0->SaveAs(output_folder + name_pdf + ".pdf");
			// cancG0->SaveAs(output_string + ".pdf");
			// cancG0->SaveAs(output_string + ".png");

			/////////////////////////////////////////
			/////////////////////////////////////////

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
			TH1D *histAccIteration = new TH1D(Form("histAccIteration%i", i), "", 50, 0.0, 2.0 * numerator_Acc / iteration_in_double);
			for (int i = 0; i < iteration_Acc; ++i)
			{
				histAccIteration->Fill(nb_jpsi_per_iteration[i]);
			}

			gStyle->SetOptStat(1110);
			histAccIteration->Draw("hist");
			histAccIteration->SetTitle(";Fitted events;Counts");

			cancAccIteration->SaveAs(output_folder + name_pdf + ".pdf");

			gStyle->SetOptStat(1);

			// Compute_Acc();

			/* //Store old code for now
			/////////////////////////////////////////
			/////Acc calculation
			/////////////////////////////////////////

			/////Simulate the background and normalize simulations to data
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

			double nb_JPsi_integral_MC_raw = h_only_JPsi->Integral(h_only_JPsi->FindBin(2.95), h_only_JPsi->FindBin(3.15));
			double normalization_MC_to_data = 2.0 * nb_JPsi / nb_JPsi_integral_MC_raw;

			cout << "Normalization factor for signal in MC: " << normalization_MC_to_data << endl;

			h_only_JPsi->Scale(normalization_MC_to_data);
			double nb_JPsi_integral_MC = h_only_JPsi->Integral();

			hs_JPsi_BG->Add(h_only_JPsi);
			hs_JPsi_BG->Add(BG_add_hist);

			// hs_JPsi->Add(BG_add_hist);
			/////End simulate the background

			TCanvas *cancAcc = new TCanvas("", "can0", 1500, 1000);
			cancAcc->cd();
			TH1D *hlast = (TH1D *)(hs_JPsi_BG->GetStack()->Last());
			hlast->SetTitle(";" + label + ";Events");
			hlast->GetYaxis()->SetTitleOffset(1.3);
			hlast->Draw("");			// hist");
			hs_JPsi_BG->Draw("e same"); // hist");

			cout << "///////////////////////" << endl;
			cout << "FIT MC" << endl;
			cout << "///////////////////////" << endl;
			Fit_Function Fit_func_MC;
			redefineErrors(hlast);
			Fit_func_MC.Set_Data_hist(hlast);
			Fit_func_MC.Set_Limits(min_fit, max_fit);
			// Fit_func_MC.Single_Gaussian_fit("SLER", Form("func_MC_%i", i));

			if (fit_procedure == "Default")
				Fit_func_MC.Single_Gaussian_Int_fit("SLER", Form("func_MC_%i", i));

			if (fit_procedure == "Crystall ball Pol 2 BG")
				Fit_func_MC.Crystall_Ball_fit("SLER", Form("func_MC_%i", i));

			if (fit_procedure == "Crystall ball exp BG")
				Fit_func_MC.Crystall_Ball_fit_exp("SLER", Form("func_MC_%i", i));

			if (fit_procedure == "Pol 2 BG")
				Fit_func_MC.Single_Gaussian_Int_fit_Pol_BG_V2("SLER", Form("func_MC_%i", i));

			if (fit_procedure == "Double Gaussian")
				Fit_func_MC.Double_Gaussian_Fit("SLR", Form("func_MC_%i", i));

			//  Fit_func_MC.Single_Gaussian_Fit_Flat_BG("SLR",Form("func_MC_%i", i));
			cout << "///////////////////////" << endl;

			// hs_JPsi_BG->SetLineWidth(1);
			// hs_JPsi_BG->SetLineColor(kBlack);
			// hs_JPsi_BG->SetMarkerSize(1);

			hlast->SetLineWidth(1);
			hlast->SetLineColor(kBlack);
			hlast->SetMarkerStyle(1);
			hlast->SetMarkerSize(1);

			hlast->Draw("");				 // hist");
			hs_JPsi_BG->Draw("e same hist"); //");//
			Fit_func_MC.Draw_Functions();

			auto legend_acc = new TLegend(0.54, 0.87, 0.90, 0.60);
			legend_acc->AddEntry(h_only_JPsi, Form("nb JPsi %3.1f", nb_JPsi_integral_MC), "f1");
			legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("J#psi fit (%3.1f #pm %3.1f)  ", Fit_func_MC.Get_Integral_Signal(), Fit_func_MC.Get_Integral_Error_Signal()), "l");
			legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("#Chi^{2} %3.1f, NdF %3.1f, #Chi^{2}/NdF %3.1f ", Fit_func_MC.chi2, Fit_func_MC.NDF, Fit_func_MC.chi2 / Fit_func_MC.NDF), "");
			legend_acc->SetFillStyle(0);
			legend_acc->SetLineWidth(0);
			legend_acc->Draw("same ");

			double nb_JPsi_MC = Fit_func_MC.Get_Integral_Signal(); // nb_JPsi_integral; //

			Acc_Num.push_back(nb_JPsi_MC / normalization_MC_to_data); //(sample_Acc->Integral());
			// cout<<"acc num 1 "<<sample_Acc->Integral()<<endl;

			cancAcc->SaveAs(output_folder + name_pdf + ".pdf");
			*/
		}
	}
	//
	// void Compute_Acc()
	//{
	//	/////////////////////////////////////////
	//	/////Acc calculation
	//	/////////////////////////////////////////
	//
	//	/////Simulate the background and normalize simulations to data
	//	float nb_event_bg = (Fit_func.function_BG->Integral(min_histo, max_histo)) / (Data_hist->GetXaxis()->GetBinWidth(2));
	//	TH1D *BG_add_hist = new TH1D("BG_add_hist", "", nbBins, min_histo, max_histo);
	//	TString name_BG_func = Form("func_%i_bg_func", i);
	//	BG_add_hist->FillRandom(name_BG_func, 50000);
	//	BG_add_hist->Scale(nb_event_bg / 50000.);
	//	// BG_add_hist->Scale(nb_event_bg/50000.); //Remove the BH background
	//	BG_add_hist->SetLineWidth(0);
	//	BG_add_hist->SetLineColor(kBlack);
	//	BG_add_hist->SetMarkerSize(0);
	//	BG_add_hist->SetFillColorAlpha(kCyan, 0.65);
	//	BG_add_hist->SetStats(kFALSE);
	//	TH1D *h_only_JPsi = (TH1D *)(hs_JPsi->GetStack()->Last());
	//
	//	double nb_JPsi_integral_MC_raw = h_only_JPsi->Integral(h_only_JPsi->FindBin(2.95), h_only_JPsi->FindBin(3.15));
	//	double normalization_MC_to_data = 2.0 * nb_JPsi / nb_JPsi_integral_MC_raw;
	//
	//	cout << "Normalization factor for signal in MC: " << normalization_MC_to_data << endl;
	//
	//	h_only_JPsi->Scale(normalization_MC_to_data);
	//	double nb_JPsi_integral_MC = h_only_JPsi->Integral();
	//
	//	hs_JPsi_BG->Add(h_only_JPsi);
	//	hs_JPsi_BG->Add(BG_add_hist);
	//
	//	// hs_JPsi->Add(BG_add_hist);
	//	/////End simulate the background
	//
	//	TCanvas *cancAcc = new TCanvas("", "can0", 1500, 1000);
	//	cancAcc->cd();
	//	TH1D *hlast = (TH1D *)(hs_JPsi_BG->GetStack()->Last());
	//	hlast->SetTitle(";" + label + ";Events");
	//	hlast->GetYaxis()->SetTitleOffset(1.3);
	//	hlast->Draw("");			// hist");
	//	hs_JPsi_BG->Draw("e same"); // hist");
	//
	//	cout << "///////////////////////" << endl;
	//	cout << "FIT MC" << endl;
	//	cout << "///////////////////////" << endl;
	//	Fit_Function Fit_func_MC;
	//	redefineErrors(hlast);
	//	Fit_func_MC.Set_Data_hist(hlast);
	//	Fit_func_MC.Set_Limits(min_fit, max_fit);
	//	// Fit_func_MC.Single_Gaussian_fit("SLER", Form("func_MC_%i", i));
	//
	//	if (fit_procedure == "Default")
	//		Fit_func_MC.Single_Gaussian_Int_fit("SLER", Form("func_MC_%i", i));
	//
	//	if (fit_procedure == "Crystall ball Pol 2 BG")
	//		Fit_func_MC.Crystall_Ball_fit("SLER", Form("func_MC_%i", i));
	//
	//	if (fit_procedure == "Crystall ball exp BG")
	//		Fit_func_MC.Crystall_Ball_fit_exp("SLER", Form("func_MC_%i", i));
	//
	//	if (fit_procedure == "Pol 2 BG")
	//		Fit_func_MC.Single_Gaussian_Int_fit_Pol_BG_V2("SLER", Form("func_MC_%i", i));
	//
	//	if (fit_procedure == "Double Gaussian")
	//		Fit_func_MC.Double_Gaussian_Fit("SLR", Form("func_MC_%i", i));
	//
	//	//  Fit_func_MC.Single_Gaussian_Fit_Flat_BG("SLR",Form("func_MC_%i", i));
	//	cout << "///////////////////////" << endl;
	//
	//	// hs_JPsi_BG->SetLineWidth(1);
	//	// hs_JPsi_BG->SetLineColor(kBlack);
	//	// hs_JPsi_BG->SetMarkerSize(1);
	//
	//	hlast->SetLineWidth(1);
	//	hlast->SetLineColor(kBlack);
	//	hlast->SetMarkerStyle(1);
	//	hlast->SetMarkerSize(1);
	//
	//	hlast->Draw("");				 // hist");
	//	hs_JPsi_BG->Draw("e same hist"); //");//
	//	Fit_func_MC.Draw_Functions();
	//
	//	auto legend_acc = new TLegend(0.54, 0.87, 0.90, 0.60);
	//	legend_acc->AddEntry(h_only_JPsi, Form("nb JPsi %3.1f", nb_JPsi_integral_MC), "f1");
	//	legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("J#psi fit (%3.1f #pm %3.1f)  ", Fit_func_MC.Get_Integral_Signal(), Fit_func_MC.Get_Integral_Error_Signal()), "l");
	//	legend_acc->AddEntry(Fit_func_MC.function_Signal, Form("#Chi^{2} %3.1f, NdF %3.1f, #Chi^{2}/NdF %3.1f ", Fit_func_MC.chi2, Fit_func_MC.NDF, Fit_func_MC.chi2 / Fit_func_MC.NDF), "");
	//	legend_acc->SetFillStyle(0);
	//	legend_acc->SetLineWidth(0);
	//	legend_acc->Draw("same ");
	//
	//	double nb_JPsi_MC = Fit_func_MC.Get_Integral_Signal(); // nb_JPsi_integral; //
	//
	//	Acc_Num.push_back(nb_JPsi_MC / normalization_MC_to_data); //(sample_Acc->Integral());
	//	// cout<<"acc num 1 "<<sample_Acc->Integral()<<endl;
	//
	//	cancAcc->SaveAs(output_folder + name_pdf + ".pdf");
	//}

	void Process_MC()
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
			double Eg_min = stof((string)labels_MC[i][7].Data());
			double Eg_max = stof((string)labels_MC[i][8].Data());

			std::vector<double> integral_flux{};
			//////////////////////////////////////////////////

			cout << "\n";
			cout << "//////////////////////////////////////////////////"
				 << "\n";
			cout << "Doing " << i << " " << label << " plot"
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
			// THStack *hs_MC_flux = new THStack("hs_MC_flux", "");
			//  THStack *hs_MC = new THStack("hs_MC", "");

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
					// TCut weight1 = Form("(%s+%s)*%f", "virtual_flux_Frixione_Gen", "real_flux_Gen", lumi_sample);
					//  TCut weight1 = Form("(%s)*%f", "virtual_flux_Gen", lumi_sample);
					TCut weight_acc = Form("%s*%f*%f*%f/(%i)", "weight", xsec, Analysis_Sample.lumi_factor, lumi_sample, nbEvents_sample); // Form("%f",1.0);//
					// TCut weight_acc = Form("%s*(%s+%s)*%f*%f*%f/((%s+%s)*%i)", "weight", "virtual_flux_Frixione_Gen", "real_flux_Gen", xsec, lumi_factor, lumi_sample, "virtual_flux_Gen", "real_flux_Gen", nbEvents_sample); // Form("%f",1.0);//

					// MC_tree[j]->Draw(label + ">>" + hist_name1, cut * weight1);
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

					// hs_MC_flux->Add(sample_hist1);

					// double integral_flux_sample = (sample_hist1->GetEntries()>100) ? (Eg_max - Eg_min) * ((sample_hist1->Integral()) / sample_hist1->GetEntries()) : 0.0;
					double integral_flux_sample = (sample_hist1->GetEntries() > 100) ? (Eg_max - Eg_min) * ((sample_hist1->Integral()) / sample_Acc->Integral()) : 0.0;
					integral_flux.push_back(integral_flux_sample);

					if (debug)
					{
						cout << " IN FLUX CALCULATION " << endl;
						cout << Analysis_Sample.samples[j][0] << "\n";
						cout << sample_hist1->Integral() << endl;
						cout << sample_hist1->GetEntries() << endl;
					}

					legend->AddEntry(sample_Acc, Form("%s (%3.6f), #sigma=%3.1f pb", Analysis_Sample.samples[j][3].Data(), sample_Acc->Integral(), stof(Analysis_Sample.samples[j][1].Data())), "f1");
				}

				if (Analysis_Sample.samples[j][3] == "Jpsi_Rad_Corr")
				{
					TString hist_no_rad = Form("sample_no_rad_%s_%i_%i", Analysis_Sample.samples[j][3].Data(), j, i);
					TH1D *sample_no_rad = new TH1D(hist_no_rad, "", nbBins, min_histo, max_histo);

					TCut weight_no_rad = Form("%s*%f*%f*%f/(%i)", "weight", xsec, Analysis_Sample.lumi_factor, lumi_sample, nbEvents_sample); // Form("%f",1.0);//

					Analysis_Sample.MC_tree[j]->Draw(label + ">>" + hist_no_rad, cut * weight_no_rad);

					// hs_MC_no_rad->SetTitle(";" + label + ";Events");
					hs_MC_no_rad->SetTitle("; M_{GEN} ;Events");
					// hs_MC_no_rad->GetYaxis()->SetTitleOffset(1.3);
					hs_MC_no_rad->Add(sample_no_rad);
				}
			}

			/////////////////////////
			// Calculation of CS
			/////////////////////////(TH1D *)(hs->GetStack()->Last())
			double avg_flux = std::accumulate(integral_flux.begin(), integral_flux.end(), 0.0);
			for (const auto &element : integral_flux)
			{
				if (debug)
					std::cout << element << " ";
			}
			cout << "Average Flux " << avg_flux << "\n";
			TH1D *MC_stack_hist = (TH1D *)(hs_MC->GetStack()->Last());
			double Acc = Acc_Num[i] / (Analysis_Sample.reduction_factor_Gen * MC_stack_hist->Integral());

			cout << "Acc calculation: num " << Acc_Num[i] << " denom: " << (MC_stack_hist->Integral()) << "\n";

			TH1D *No_rad_MC_stack_hist = (TH1D *)(hs_MC_no_rad->GetStack()->Last());
			double Rad_corr = (MC_stack_hist->Integral()) / (No_rad_MC_stack_hist->Integral());

			Acc_hist->SetBinContent(i + 1, Acc);
			Flux_hist->SetBinContent(i + 1, avg_flux * Analysis_Sample.lumi_factor);
			Nb_JPsi_hist->SetBinContent(i + 1, nb_JPsi_Data[i]);
			Rad_corr_hist->SetBinContent(i + 1, Rad_corr);

			cout << "Acc " << Acc << "\n";
			cout << "Rad_corr " << Rad_corr << "\n";
			cout << "Nb JPsi " << nb_JPsi_Data[i] << "\n";
			cout << "Br " << Branching_ratio << "\n";
			cout << "w_c " << w_c[i] << "\n";
			cout << "Size of the bin Delta_E " << (Eg_max - Eg_min) << "\n";

			double CS_usual = 0.001 * nb_JPsi_Data[i] / (avg_flux * Analysis_Sample.lumi_factor * Branching_ratio * Rad_corr * (Acc)*w_c[i]);
			cout << "CS using normal formula " << CS_usual << "\n"; // remove lumi factor ?

			double error_CS_usual = 0.001 * error_nb_JPsi_Data[i] / (avg_flux * Analysis_Sample.lumi_factor * Branching_ratio * Rad_corr * (Acc)*w_c[i]); //
			cout << "error on the jpsi number " << error_nb_JPsi_Data[i] << "\n";
			cout << "error CS using normal formula " << error_CS_usual << "\n";

			cout << "average variable " << average_variable[i] << "\n";

			Latex_Table.Add_value(i, labels_MC.size(), CS_usual, error_CS_usual, average_variable[i], sigma_variable[i]);

			JPsi_CS_Graph.SetPoint(i, average_variable[i], CS_usual);
			// JPsi_CS_Graph.SetPointError(i, average_variable[i] - Eg_min, Eg_max - average_variable[i], error_CS_usual, error_CS_usual);
			JPsi_CS_Graph.SetPointError(i, sigma_variable[i], sigma_variable[i], error_CS_usual, error_CS_usual);

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

			gSystem->ProcessEvents();
		}

		// Analysis_Sample.clearAttribute();
		cout << "Size " << Analysis_Sample.MC_tree.size() << endl;
	}

	void Finalize_int_CS_Calculation()
	{
		//////////////////// Glue X and Hall C results ////////////////////////////////
		TGraphErrors *JPsi_CS_GlueX = new TGraphErrors(18, EgGluex, xsGluex, xeGluex, sxsGluex);
		JPsi_CS_GlueX->SetMarkerColor(kRed);
		JPsi_CS_GlueX->SetLineColor(kRed);
		JPsi_CS_GlueX->SetMarkerStyle(22);
		JPsi_CS_GlueX->SetMarkerSize(3.);

		TGraphErrors *JPsi_CS_GlueX_old = new TGraphErrors(10, EgGluex_old, xsGluex_old, xeGluex_old, sxsGluex_old);
		JPsi_CS_GlueX_old->SetMarkerColor(kRed + 4);
		JPsi_CS_GlueX_old->SetLineColor(kRed + 4);
		JPsi_CS_GlueX_old->SetMarkerStyle(22);
		JPsi_CS_GlueX_old->SetMarkerSize(3.);

		TGraphErrors *JPsi_CS_HallC = new TGraphErrors(10, Ei, xstot, xetot, sxstot);
		JPsi_CS_HallC->SetMarkerColor(kGreen);
		JPsi_CS_HallC->SetLineColor(kGreen);
		JPsi_CS_HallC->SetMarkerStyle(23);
		JPsi_CS_HallC->SetMarkerSize(3.);

		TGraphErrors *Richard_graph = new TGraphErrors(8, E_Richard, xstot_Richard, s_E_Richard, s_xstot_Richard);
		Richard_graph->SetMarkerColor(kGreen - 2);
		Richard_graph->SetLineColor(kGreen - 2);
		Richard_graph->SetMarkerStyle(23);
		Richard_graph->SetMarkerSize(3.);

		TF1 *f_JPsi_sigm_int_2g = new TF1("f_JPsi_sigm_int_2g", JPsi_sigm_int_2g, 8.25, 10.7, 3);
		double SLAC_Fit_scale = 7.79117e-23;
		double tSlope = 1.13;
		f_JPsi_sigm_int_2g->SetParameter(1, SLAC_Fit_scale);
		f_JPsi_sigm_int_2g->SetParameter(2, tSlope);

		/*TFile *CS_file = new TFile("../Spring2019/outputTCS_JPsi_Spring2019_Pass2.root");
		TTree *CS_tree = (TTree *)CS_file->Get("tree_Gen");
		TH1F *CS_hist_1 = new TH1F("CS_hist_1", "CS_hist_1", 50, 8.6, 10.2);
		CS_tree->Draw("Epho_Gen>>CS_hist_1", "weight*(1.0)/(((10.2-8.6)/50.)*1124999*(virtual_flux_Gen+real_flux_Gen))", "hist");*/
		// TFile *CS_file = new TFile(Analysis_Sample.sample_folder + "Simulation/JPsi_Rad_corr_Fall2018_45_022024.root");
		// TFile *CS_file = new TFile(Analysis_Sample.sample_folder + "JPsi_RGB_Spring19.root");
		// TTree *CS_tree = (TTree *)CS_file->Get("tree_Gen");
		// TH1F *CS_hist_1 = new TH1F("CS_hist_1", "CS_hist_1", 50, 8.6, 10.6);
		// CS_tree->Draw("Epho_Gen>>CS_hist_1", "weight*(1.0)/(((10.6-8.6)/50.)*2015000*(virtual_flux_Gen+real_flux_Gen))", "hist");

		/*auto CS_graph_theo = new TGraph();
		for (int i = 1; i <= CS_hist_1->GetNbinsX(); i++)
		{
			CS_graph_theo->SetPoint(i - 1, CS_hist_1->GetBinCenter(i), CS_hist_1->GetBinContent(i));
			cout << CS_hist_1->GetBinCenter(i) << " " << CS_hist_1->GetBinContent(i) << " " << CS_hist_1->GetNbinsX() << endl;
		}*/
		///////////////////////////////////////////////////////////////////////////////

		TCanvas *cancG0 = new TCanvas("", "can0", 1500, 1000);

		JPsi_CS_GlueX->SetTitle(";E_{#gamma} [GeV]; #sigma [nb]");
		JPsi_CS_GlueX->GetXaxis()->SetRangeUser(8., 11.5);

		// JPsi_CS_Graph.SetMaximum(1.);
		JPsi_CS_GlueX->SetMinimum(0.);
		JPsi_CS_GlueX->SetMaximum(2.1);
		JPsi_CS_GlueX->Draw("AP");
		// JPsi_CS_HallC->Draw("P");
		// JPsi_CS_GlueX_old->Draw("P");
		if (RunGroup == "RGB")
			Richard_graph->Draw("P");
		JPsi_CS_Graph.Draw("P");
		// f_JPsi_sigm_int_2g->Draw("same");
		// CS_graph_theo->Draw("same");
		auto legend_CS = new TLegend(0.12, 0.7, 0.48, 0.85);
		legend_CS->AddEntry(&JPsi_CS_Graph, "CLAS12 (this study)", "lp");
		// legend_CS->AddEntry(JPsi_CS_HallC, "Hall C (unpublished)", "lp");
		legend_CS->AddEntry(JPsi_CS_GlueX, "GlueX (2023)", "lp");
		// legend_CS->AddEntry(JPsi_CS_GlueX_old, "GlueX (2019)", "lp");
		if (RunGroup == "RGB")
			legend_CS->AddEntry(Richard_graph, "Richard's analysis", "lp");
		legend_CS->SetFillStyle(0);
		legend_CS->SetLineWidth(0);
		legend_CS->Draw("same");

		TText *t = new TText();
		t->SetTextAlign(22);
		// t->SetTextColorAlpha(kGray, 0.70);
		t->SetTextColor(kGray);
		t->SetTextFont(40);
		t->SetTextSize(0.25);
		t->SetTextAngle(25);
		t->DrawTextNDC(.5, .53, "Preliminary");

		cancG0->SaveAs(output_folder + name_pdf + ".pdf");

		/// Log CS
		TCanvas *cancG01 = new TCanvas("", "can01", 1500, 1000);

		gPad->SetLogy();
		JPsi_CS_GlueX->SetMaximum(10.0);
		JPsi_CS_GlueX->Draw("AP");
		// JPsi_CS_HallC->Draw("P");
		// JPsi_CS_GlueX_old->Draw("P");
		JPsi_CS_Graph.Draw("P");
		// f_JPsi_sigm_int_2g->Draw("same");
		// CS_graph_theo->Draw("same");

		legend_CS->Draw("same");

		t->DrawTextNDC(.5, .53, "Preliminary");

		cancG01->SaveAs(output_folder + name_pdf + ".pdf");

		/// Blind CS
		TCanvas *cancG1 = new TCanvas("", "can1", 1500, 1000);

		// JPsi_CS_Graph.SetMaximum(1.);
		auto JPsi_CS_Graph_blind = new TGraphErrors(6);
		for (int i = 0; i < labels_MC.size(); i++)
		{
			JPsi_CS_Graph_blind->SetPoint(i, JPsi_CS_Graph.GetPointX(i), 0.8);
			JPsi_CS_Graph_blind->SetPointError(i, 0.0, (JPsi_CS_Graph.GetErrorY(i)));
		}
		JPsi_CS_Graph_blind->SetTitle(";E_{#gamma} [GeV];#sigma [nb]");
		JPsi_CS_GlueX->SetTitle(";E_{#gamma} [GeV];#sigma [nb]");
		JPsi_CS_GlueX->GetXaxis()->SetRangeUser(8., 11.5);
		JPsi_CS_Graph_blind->SetMarkerColor(kBlue);
		JPsi_CS_Graph_blind->SetLineColor(kBlue);
		JPsi_CS_Graph_blind->SetMarkerStyle(21);
		JPsi_CS_Graph_blind->SetMarkerSize(2.);
		JPsi_CS_Graph_blind->SetMinimum(0.);
		JPsi_CS_Graph_blind->SetMaximum(1.8);
		JPsi_CS_HallC->SetMinimum(0.);
		JPsi_CS_HallC->SetMaximum(1.8);
		// JPsi_CS_HallC->Draw("P");
		JPsi_CS_GlueX->SetMaximum(3.0);
		JPsi_CS_GlueX->Draw("AP");
		// JPsi_CS_GlueX_old->Draw("P");
		JPsi_CS_Graph_blind->Draw("P");
		// f_JPsi_sigm_int_2g->Draw("same");
		// CS_graph_theo->Draw("same");
		auto legend_CS_1 = new TLegend(0.12, 0.7, 0.48, 0.85);
		legend_CS_1->AddEntry(JPsi_CS_Graph_blind, "Projected stat. error bars", "lp");
		// legend_CS->AddEntry(JPsi_CS_HallC, "Hall C (unpublished)", "lp");
		legend_CS_1->AddEntry(JPsi_CS_GlueX, "GlueX (2023)", "lp");
		// legend_CS_1->AddEntry(JPsi_CS_GlueX_old, "GlueX (2019)", "lp");
		legend_CS_1->SetFillStyle(0);
		legend_CS_1->SetLineWidth(0);
		legend_CS_1->Draw("same");

		t->DrawTextNDC(.5, .53, "Preliminary");

		cancG1->SaveAs(output_folder + name_pdf + ".pdf");

		gStyle->SetPaintTextFormat("4.3f");

		TCanvas *cancG2 = new TCanvas("", "can2", 1500, 1000);
		Acc_hist->Draw("hist text0");
		cancG2->SaveAs(output_folder + name_pdf + ".pdf");

		TCanvas *cancG3 = new TCanvas("", "can3", 1500, 1000);
		Flux_hist->Draw("hist text0");
		cancG3->SaveAs(output_folder + name_pdf + ".pdf");

		TCanvas *cancG5 = new TCanvas("", "can5", 1500, 1000);
		Nb_JPsi_hist->Draw("hist text0");
		cancG5->SaveAs(output_folder + name_pdf + ".pdf");

		TCanvas *cancG4 = new TCanvas("", "can4", 1500, 1000);
		Rad_corr_hist->Draw("hist text0");
		cancG4->SaveAs(output_folder + name_pdf + ".pdf)");
	}
	//////////////////////////////

	/////// Processing function for the differential CS //////
	void Process_REC_differential()
	{
	}

	void Process_MC_differential()
	{
	}
	//////////////////////////////
};

#endif
