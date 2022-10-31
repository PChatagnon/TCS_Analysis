#ifndef TCSParameters
#define TCSParameters

#include <iostream>
#include <fstream>

class Parameters
{
public:
        float MMassBeamCut;
        float PtCut;
        float ElecMomCut;
        float PosiMomCut;
        float ElecMinSF;
        float PosiMinSF;
        float MLPscoreCut;
        float MLPMomCut;
        bool RadCorr;
        bool CentralMomCorr;
        float factorPola;

        Parameters(string inputfile)
        {
                std::ifstream infile(inputfile);

                string value;
                string cut;

                while (infile >> cut >> value)
                {
                        if (cut == "MMassBeamCut")
                                MMassBeamCut = stof(value);
                        if (cut == "PtCut")
                                PtCut = stof(value);
                        if (cut == "ElecMomCut")
                                ElecMomCut = stof(value);
                        if (cut == "PosiMomCut")
                                PosiMomCut = stof(value);
                        if (cut == "ElecMinSF")
                                ElecMinSF = stof(value);
                        if (cut == "PosiMinSF")
                                PosiMinSF = stof(value);
                        if (cut == "MLPscoreCut")
                                MLPscoreCut = stof(value);
                        if (cut == "MLPMomCut")
                                MLPMomCut = stof(value);
                        if (cut == "RadCorr")
                                RadCorr = (value == "true");
                        if (cut == "CentralMomCorr")
                                CentralMomCorr = (value == "true");
                        if (cut == "factorPola")
                                factorPola = stof(value);
                }
        }
};

#endif
