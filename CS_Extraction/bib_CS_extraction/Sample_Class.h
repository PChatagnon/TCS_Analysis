#ifndef Sample_Class
#define Sample_Class

class Sample
{
public:
        // By default RGA values are used
        TString sample_folder = "/Users/pc281174/Desktop/JPsi_analysis/Pass2_Samples/";
        float lumi_factor = 1316.875;

        // Store the samples
        std::vector<vector<TString>> samples{};

        std::vector<TTree *> reduced_samples_tree{};
        std::vector<int> ngen{};
        std::vector<TTree *> MC_tree{};
        double reduction_factor_Gen = 2.0;

        TChain *Data_tree = new TChain("tree");
        TTree *filtered_Data_tree;

        // Store color of the histogramms
        TString color_JPsi = "420";
        TString color_TCSGen = "410";
        TString color_Grape = "412";
        TString color_BG = "42";


        TString charge_45_fall_18_inbending = "26.312";
        TString charge_50_fall_18_inbending = "4.0006";
        TString charge_55_fall_18_inbending = "5.35578";

        TString charge_40_fall_18_outbending = "11.8306";
        TString charge_50_fall_18_outbending = "20.6199";

        TString charge_50_spring_19_inbending = "45.993";


        Sample() {}

        void clearAttribute()
        {
                //delete Data_tree;
                //delete filtered_Data_tree;
//
                //for (TTree *tree : reduced_samples_tree)
                //{
                //        delete tree; // Release memory for each TTree
                //}
                //reduced_samples_tree.clear();
//
                //for (TTree *tree : MC_tree)
                //{
                //        delete tree; // Release memory for each TTree
                //}
                //MC_tree.clear();
        }

        void Setup_RGA(TString mom_corr_version)
        {
                float Norm_factor = 2.00; // 2.00; // This factor cancels out everywhere. Could be removed
                float cs_Jpsi = 57.6;
                TString string_cs_Jpsi = to_string(Norm_factor * cs_Jpsi);
                TString string_cs_Jpsi_test = to_string(Norm_factor * cs_Jpsi * 1.5);
                TString cs_no_rad_2018 = to_string(Norm_factor * 138.24); // Normalization to 57.6*Eg_psf(10.6-8.2)
                TString cs_no_rad_2019 = to_string(Norm_factor * 115.2);  // Normalization to 57.6*Eg_psf(10.2-8.2)

               
                // OLD Charges
                /*TString charge_45_fall_18_inbending = "34.1137765";
                TString charge_50_fall_18_inbending = "3.4429876";
                TString charge_55_fall_18_inbending = "5.8921636";

                TString charge_40_fall_18_outbending = "12.013926";
                TString charge_50_fall_18_outbending = "21.8007092";

                TString charge_50_spring_19_inbending = "50.5319";*/

                // JPsi BH radiative correction
                //  Fall2018
                //        samples.push_back({sample_folder + "Simulation/JPsi_Rad_corr_Fall2018_45_022024.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_45_fall_18_inbending});
                //        samples.push_back({sample_folder + "Simulation/JPsi_Rad_corr_Fall2018_55_022024.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_55_fall_18_inbending});
                //        samples.push_back({sample_folder + "Simulation/JPsi_Rad_corr_Fall2018_50_022024.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_inbending});
                //        // Outbending
                //        samples.push_back({sample_folder + "Simulation/JPsi_Rad_corr_Fall2018_40_022024.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_40_fall_18_outbending});
                //        samples.push_back({sample_folder + "Simulation/JPsi_Rad_corr_Fall2018_50_out_022024.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_outbending});
                //        // Spring2019
                //        samples.push_back({sample_folder + "Simulation/JPsi_Rad_corr_Spring2019_022024.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});



                //////Review version 1//////
/*
                // JPsi VM radiative correction
                
                 //  Fall2018
                samples.push_back({sample_folder + "Simu_Mom_Corr/JPsi_VM_Fall2018_45.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_45_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_Mom_Corr/JPsi_VM_Fall2018_55.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_55_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_Mom_Corr/JPsi_VM_Fall2018_50.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simu_Mom_Corr/JPsi_VM_Fall2018_40_out.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simu_Mom_Corr/JPsi_VM_Fall2018_50_out.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_outbending});
                // Spring2019
                samples.push_back({sample_folder + "Simu_Mom_Corr/JPsi_VM_Spring2019.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});
                //samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_Spring2019_SDME_s.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});



                ////////Review version 1//////
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_45_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_45_fall_18_inbending}); // Normalization to 4*57.6*Eg_psf(10.6-8.2)
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_50_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_inbending});
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_55_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_55_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_40_out_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_50_out_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_outbending});
                // Spring2019
                samples.push_back({sample_folder + "Simulation/JPsi_Spring2019_022024.root", cs_no_rad_2019, color_JPsi, "Jpsi_Rad_Corr", charge_50_spring_19_inbending});
                
*/



                // JPsi VM radiative correction + SDME //round 1 comments
                 //  Fall2018
                /*samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_45nA_SDME.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_45_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_55nA_SDME.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_55_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_50nA_SDME.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_40nA_out_SDME.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_50nA_out_SDME.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_outbending});
                // Spring2019
                //samples.push_back({sample_folder + "Simu_Mom_Corr/JPsi_VM_Spring2019.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});
                samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_Spring2019_new_June25.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});
                //samples.push_back({sample_folder + "Simu_JPsi_SDME/Frixione_25/outputTCS_Spring2019_SDME.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});*/

                /* //////No rad + SDME//////
                samples.push_back({sample_folder + "Simu_JPsi_SDME/outputTCS_45nA_SDME_NoRad.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_45_fall_18_inbending}); // Normalization to 4*57.6*Eg_psf(10.6-8.2)
                samples.push_back({sample_folder + "Simu_JPsi_SDME/outputTCS_50nA_SDME_NoRad.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_JPsi_SDME/outputTCS_55nA_SDME_NoRad.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_55_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simu_JPsi_SDME/outputTCS_40nA_out_SDME_NoRad.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simu_JPsi_SDME/outputTCS_50nA_out_SDME_NoRad.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_outbending});
                // Spring2019
                samples.push_back({sample_folder + "Simu_JPsi_SDME/outputTCS_Spring2019_SDME_NoRad.root", cs_no_rad_2019, color_JPsi, "Jpsi_Rad_Corr", charge_50_spring_19_inbending});*/
                // JPsi VM radiative correction + SDME //round 1 comments


                // JPsi VM radiative correction + SDME + variation flux //September 2025
                //  Fall2018
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_JPsi_Flux_Fall2018_45.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_45_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_JPsi_Flux_Fall2018_55.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_55_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_JPsi_Flux_Fall2018_50.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_JPsi_Flux_Fall2018_40_out.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_JPsi_Flux_Fall2018_50_out.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_outbending});
                // Spring2019
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_JPsi_Flux_Spring2019.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});


                //////No rad + SDME + Flux//////
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_Flux_No_Rad_45.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_45_fall_18_inbending}); // Normalization to 4*57.6*Eg_psf(10.6-8.2)
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_Flux_No_Rad_50.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_inbending});
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_Flux_No_Rad_55.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_55_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_Flux_No_Rad_40_out.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_Flux_No_Rad_50_out.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_outbending});
                // Spring2019
                samples.push_back({sample_folder + "Simu_JPsi_Flux/outputTCS_Flux_No_Rad_Spring19.root", cs_no_rad_2019, color_JPsi, "Jpsi_Rad_Corr", charge_50_spring_19_inbending});
                // JPsi VM radiative correction + SDME + variation flux //September 2025


                
                //samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_45_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_45_fall_18_inbending}); // Normalization to 4*57.6*Eg_psf(10.6-8.2)
                //samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_50_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_inbending});
                //samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_55_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_55_fall_18_inbending});
                //// Outbending
                //samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_40_out_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_40_fall_18_outbending});
                //samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_50_out_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_outbending});
                //// Spring2019
                //samples.push_back({sample_folder + "Simulation/JPsi_Spring2019_022024.root", cs_no_rad_2019, color_JPsi, "Jpsi_Rad_Corr", charge_50_spring_19_inbending});






/*
                //No mom. Corr
                //  Fall2018
                samples.push_back({sample_folder + "Simulation/JPsi_VM_Rad/JPsi_VM_Fall2018_45.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_45_fall_18_inbending});
                samples.push_back({sample_folder + "Simulation/JPsi_VM_Rad/JPsi_VM_Fall2018_55.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_55_fall_18_inbending});
                samples.push_back({sample_folder + "Simulation/JPsi_VM_Rad/JPsi_VM_Fall2018_50.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simulation/JPsi_VM_Rad/JPsi_VM_Fall2018_40_out.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simulation/JPsi_VM_Rad/JPsi_VM_Fall2018_50_out.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_fall_18_outbending});
                // Spring2019
                samples.push_back({sample_folder + "Simulation/JPsi_VM_Rad/JPsi_VM_Spring2019.root", string_cs_Jpsi, color_JPsi, "J#psi", charge_50_spring_19_inbending});
                


                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_45_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_45_fall_18_inbending});
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_50_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_inbending});
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_55_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_55_fall_18_inbending});
                // Outbending
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_40_out_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_40_fall_18_outbending});
                samples.push_back({sample_folder + "Simulation/JPsi_Fall2018_50_out_022024.root", cs_no_rad_2018, color_JPsi, "Jpsi_Rad_Corr", charge_50_fall_18_outbending});
                // Spring2019
                samples.push_back({sample_folder + "Simulation/JPsi_Spring2019_022024.root", cs_no_rad_2019, color_JPsi, "Jpsi_Rad_Corr", charge_50_spring_19_inbending});
*/
                
                /*
                                          // TCSGEN
                                          //  Fall2018
                                          samples.push_back({sample_folder + "Simulation/TCS_Gen_Rad_corr_Fall2018_45_022024.root", "1", color_TCSGen, "BH TCSGen", charge_45_fall_18_inbending});
                                          samples.push_back({sample_folder + "Simulation/TCS_Gen_Rad_corr_Fall2018_50_022024.root", "1", color_TCSGen, "BH TCSGen", charge_50_fall_18_inbending});
                                          samples.push_back({sample_folder + "Simulation/TCS_Gen_Rad_corr_Fall2018_55_022024.root", "1", color_TCSGen, "BH TCSGen", charge_55_fall_18_inbending});
                                          // Outbending
                                          samples.push_back({sample_folder + "Simulation/TCS_Gen_Rad_corr_Fall2018_40_022024.root", "1", color_TCSGen, "BH TCSGen", charge_40_fall_18_outbending});
                                          samples.push_back({sample_folder + "Simulation/TCS_Gen_Rad_corr_Fall2018_50_out_022024.root", "1", color_TCSGen, "BH TCSGen", charge_50_fall_18_outbending});
                                          // Spring2019
                                          samples.push_back({sample_folder + "Simulation/TCS_Gen_Rad_corr_Spring2019_022024.root", "1", color_TCSGen, "BH TCSGen", charge_50_spring_19_inbending});

                                          // Grape
                                          //  Fall2018
                                          samples.push_back({sample_folder + "Simulation/Grape_Rad_corr_Fall2018_45_022024.root", "0.412", color_Grape, "BH Grape", charge_45_fall_18_inbending});
                                          samples.push_back({sample_folder + "Simulation/Grape_Rad_corr_Fall2018_50_022024.root", "0.412", color_Grape, "BH Grape", charge_50_fall_18_inbending});
                                          samples.push_back({sample_folder + "Simulation/Grape_Rad_corr_Fall2018_55_022024.root", "0.412", color_Grape, "BH Grape", charge_55_fall_18_inbending});
                                          // Outbending
                                          samples.push_back({sample_folder + "Simulation/Grape_Rad_corr_Fall2018_40_022024.root", "0.412", color_Grape, "BH Grape", charge_40_fall_18_outbending});
                                          samples.push_back({sample_folder + "Simulation/Grape_Rad_corr_Fall2018_50_out_022024.root", "0.412", color_Grape, "BH Grape", charge_50_fall_18_outbending});
                                          // Spring2019
                                          samples.push_back({sample_folder + "Simulation/Grape_Rad_corr_Spring2019_022024.root", "0.320", color_Grape, "BH Grape", charge_50_spring_19_inbending});
                          */

                /*TString data_file_adress_0 = sample_folder + "Data/Data_pass2_spring2019_inbending.root";
                TString data_file_adress_1 = sample_folder + "Data/Data_pass2_fall2018_inbending.root";
                TString data_file_adress_2 = sample_folder + "Data/Data_pass2_fall2018_outbending.root";*/

                TString data_file_adress_0;
                TString data_file_adress_1;
                TString data_file_adress_2;

                if(mom_corr_version=="Default"){
                        data_file_adress_0 = sample_folder + "Data_Mom_Corr/Data_pass2_spring2019_inbending_mom_corr.root";
                        data_file_adress_1 = sample_folder + "Data_Mom_Corr/Data_pass2_fall2018_inbending_mom_corr.root";
                        data_file_adress_2 = sample_folder + "Data_Mom_Corr/Data_pass2_fall2018_outbending_mom_corr.root";
                }
                else if(mom_corr_version=="2_1"){
                        cout << "Using mom_corr version 2.1" << endl;
                        data_file_adress_0 = "/Users/pc281174/Desktop/JPsi_analysis/Mom_corr_syst/Data_pass2_fall2018_inbending_momcorr_v21.root";
                        data_file_adress_1 = "/Users/pc281174/Desktop/JPsi_analysis/Mom_corr_syst/Data_pass2_spring2019_inbending_momcorr_v21.root";
                        data_file_adress_2 = "/Users/pc281174/Desktop/JPsi_analysis/Mom_corr_syst/Data_pass2_fall2018_outbending_momcorr_v21.root";
                
                }
                else if(mom_corr_version=="2_2"){
                        data_file_adress_0 = "/Users/pc281174/Desktop/JPsi_analysis/Mom_corr_syst/Data_pass2_fall2018_inbending_momcorr_v22.root";
                        data_file_adress_1 = "/Users/pc281174/Desktop/JPsi_analysis/Mom_corr_syst/Data_pass2_spring2019_inbending_momcorr_v22.root";
                        data_file_adress_2 = "/Users/pc281174/Desktop/JPsi_analysis/Mom_corr_syst/Data_pass2_fall2018_outbending_momcorr_v22.root";
                }

                // Event mixing Background
                // samples.push_back({"/mnt/c/Users/pierrec/Desktop/ABCD_Method/BG_inbending_Spring2019_Q2_20_Formatted_3D_Fall2018.root", "0.050", color_BG, "BG", "0.00075"});
                // samples.push_back({"/mnt/c/Users/pierrec/Desktop/ABCD_Method/BG_inbending_4_Q05_Epho8_M2_Formatted_2_3D.root", "0.050", color_BG, "BG", "0.00075"});
                // samples.push_back({"/mnt/c/Users/pierrec/Desktop/ABCD_Method/BG_outbending_Q2_20_Formatted_3D_Fall2018.root", "0.050", color_BG, "BG", "0.00075"});

                cout << "Add first data file" << endl;
                Data_tree->Add(data_file_adress_0);
                cout << "Add second data file" << endl;
                Data_tree->Add(data_file_adress_1);
                cout << "Add third data file" << endl;
                Data_tree->Add(data_file_adress_2);
        }

        void Setup_RGB()
        {
                // JPsi
                //  Fall2018
                sample_folder = "/mnt/c/Users/pierrec/Desktop/TCS_Analysis/TCS_Analysis_2022/TCS_Analysis/Pass2_Samples_RGB/";

                lumi_factor = 1514.34;

                samples.push_back({sample_folder + "JPsi_RGB_Spring19.root", "0.07", color_JPsi, "J#psi", "27.237"});
                samples.push_back({sample_folder + "JPsi_RGB_Spring19_10_2.root", "0.07", color_JPsi, "J#psi", "39.389"});

                samples.push_back({sample_folder + "JPsi_RGB_Spring19.root", "0.07", color_JPsi, "Jpsi_Rad_Corr", "27.237"});
                samples.push_back({sample_folder + "JPsi_RGB_Spring19_10_2.root", "0.07", color_JPsi, "Jpsi_Rad_Corr", "39.389"});

                TString data_file_adress_0 = sample_folder + "Data_pass2_RGB_Spring19.root";
                TString data_file_adress_1 = sample_folder + "Data_pass2_RGB_Spring19_10_2.root";

                cout << "Add first data file" << endl;
                Data_tree->Add(data_file_adress_0);
                cout << "Add second data file" << endl;
                Data_tree->Add(data_file_adress_1);
        }

        void Setup_RGA_inbending_2018()
	{
		// JPsi
		//  Fall2018 Inbending

                TString data_file_adress_0 = sample_folder + "Data_Mom_Corr/Data_pass2_spring2019_inbending_mom_corr.root";
                TString data_file_adress_1 = sample_folder + "Data_Mom_Corr/Data_pass2_fall2018_inbending_mom_corr.root";
                TString data_file_adress_2 = sample_folder + "Data_Mom_Corr/Data_pass2_fall2018_outbending_mom_corr.root";

		Data_tree->Add(data_file_adress_0);
                Data_tree->Add(data_file_adress_1);
                Data_tree->Add(data_file_adress_2);
		
	}

        void Filter_Sample_Trees(TCut kinematic_cut, TCut kinematic_cut_BG, TCut exclusivity_cut)
        {
                TRandom3 random(static_cast<unsigned int>(std::time(0)));
                cout << "Start reducing sample trees... \n";
                for (int j = 0; j < samples.size(); j++)
                {
                        cout << "File " << (j + 1) << " out of " << samples.size() << endl;
                        TFile *sample_file = new TFile(samples[j][0]);
                        int nbEvents_sample = 0;
                        if (samples[j][3] == "BG")
                        {
                                nbEvents_sample = 1;
                        }
                        else
                                nbEvents_sample = ((TH1D *)sample_file->Get("evt_count"))->GetBinContent(1);
                        ngen.push_back(nbEvents_sample);

                        TTree *sample_tree = (TTree *)sample_file->Get("tree");
                        TFile *f2 = new TFile("small.root", "recreate");
                        TTree *filtered_sample_tree = sample_tree->CopyTree(kinematic_cut * exclusivity_cut);
                        if (samples[j][3] == "BG")
                        {
                                filtered_sample_tree = sample_tree->CopyTree(kinematic_cut_BG * exclusivity_cut);
                        }
                        reduced_samples_tree.push_back(filtered_sample_tree);

                        TTree *MC_sample_tree = (TTree *)sample_file->Get("tree_Gen");
                        TFile *f3 = new TFile("small2.root", "recreate");
                        TString selection_MC = Form("((virtual_flux_Gen)>0) && rndm() < %f", (1.0 / reduction_factor_Gen));
                        ; //

                        TTree *filtered_MC_tree = MC_sample_tree->CopyTree(selection_MC); //(TTree *)sample_file->Get("tree_Gen");

                        if (samples[j][3] == "J#psi")
                        {
                                filtered_MC_tree = MC_sample_tree->CopyTree(selection_MC); //"((virtual_flux_Gen)>0)");
                        }
                        cout << filtered_MC_tree->GetEntries() << endl;
                        MC_tree.push_back(filtered_MC_tree);
                }
                cout << "Finished reducing sample trees !\n";
        }

        void Filter_Data_Tree(TCut kinematic_cut, TCut data_cut, TCut exclusivity_cut)
        {
                cout << "Filtering data now...\n";

                TFile *Data_file_temp = new TFile("smalldata.root", "recreate");
                filtered_Data_tree = (TTree *)Data_tree->CopyTree(kinematic_cut * exclusivity_cut * data_cut);

                cout << "Finished reducing data !\n";
        }
};

#endif
