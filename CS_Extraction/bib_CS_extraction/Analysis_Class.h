#ifndef Analysis_Class
#define Analysis_Class

class Analysis
{
public:

	//Store int and double variables
	float Branching_ratio = 0.06;
	float lumi_factor = 1316.875;

	//Store configuration
	TString RunGroup = "RGA";
	TString variable = "Epho";
	bool debug = false;
	bool ratio_pad = false;

	//Store fit and plotting configuration
	double min_fit = 2.7;
	double max_fit = 3.3;

	//Store normalization configuration
	double Mass_norm_low = 2.6;
	double Mass_norm_high = 2.9;

	
	//Store all quantities for the analysis
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
	
	//Store the path to data and output
	TString folder_pass2 = "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/Pass2_Samples/";
	string latex_output_folder = "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/CS_Extraction";
	TString name_pdf = "CS_Extraction_CrystalBall_exp_1";

	//Store kinematic and exclusivity cuts
	TCut kinematic_cut = "pass_EC_cut &&  Proton.Theta()*180./3.141592<35. && M>2.6 && (Electron.P() > 1.7) && (Positron.P() > 1.7) && positron_SF>0.15 && electron_SF>0.15 && ( Positron.P()<4.0 || (Positron.P()>4.0 && positron_score>0.05)) && ( Electron.P()<4.0 || (Electron.P()>4.0 && electron_score>0.05))  && positron_HTCC_ECAL_match==1. && electron_HTCC_ECAL_match==1.";
	TCut kinematic_cut_BG = "Proton.Theta()*180./3.141592<35. &&  M>2.6 && (Electron.P() > 1.7) && (Positron.P() > 1.7)";
	TCut exclusivity_cut = "abs(MMassBeam)<0.4  && abs(Q2)<0.5";
	TCut data_cut = "abs(positron_HTCCt-electron_HTCCt)<4";

	Analysis() {}

};

#endif
