#!/bin/bash

#Standard
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Standard_config_bin1.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Standard_config_bin2.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Standard_config_bin3.dat
#
##Q2 cut
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_02_config_bin1.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_02_config_bin2.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_02_config_bin3.dat
#
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_08_config_bin1.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_08_config_bin2.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_08_config_bin3.dat
#
#
##MM cut
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/MM_02_config_bin1.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/MM_02_config_bin2.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/MM_02_config_bin3.dat

#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/MM_08_config_bin1.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/MM_08_config_bin2.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/MM_08_config_bin3.dat


#Fit procedure

#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_02_config_bin1.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_02_config_bin2.dat
root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_02_config_bin3.dat
#
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_08_config_bin1.dat
#root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_08_config_bin2.dat
root -l CS_Extraction_combine_Diff_t_cs.C -param Results_diff_CS/Q2_08_config_bin3.dat