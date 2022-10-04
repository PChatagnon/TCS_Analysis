#ifndef TCSPositronClass
#define TCSPositronClass

class PositronIdentification
{
public:
        TString methodName;
        TString weightFile;
        double scoreCut;
        double momentumThreshold;
        Float_t SFPCAL, SFECIN, SFECOUT;
        Float_t m2PCAL, m2ECIN, m2ECOUT;
        TMVA::Reader *PositronIdentificationReader;
        double score;

        PositronIdentification(TString input_methodName, TString input_weightFile, double input_scoreCut, double input_momentumThreshold)
        {
                methodName = input_methodName;
                weightFile = input_weightFile;
                PositronIdentificationReader = new TMVA::Reader("!Color:!Silent");
                scoreCut = input_scoreCut;
                momentumThreshold = input_momentumThreshold;

                score = -1.;
        }

        void AddVariable(Float_t input_variable, TString input_variable_string)
        {
                PositronIdentificationReader->AddVariable(input_variable_string, &input_variable);
        }

        void InitializePositronIdentification()
        {
                PositronIdentificationReader->AddVariable("SFPCAL", &SFPCAL);
                PositronIdentificationReader->AddVariable("SFECIN", &SFECIN);
                PositronIdentificationReader->AddVariable("SFECOUT", &SFECOUT);
                PositronIdentificationReader->AddVariable("m2PCAL", &m2PCAL);
                PositronIdentificationReader->AddVariable("m2ECIN", &m2ECIN);
                PositronIdentificationReader->AddVariable("m2ECOUT", &m2ECOUT);
                PositronIdentificationReader->BookMVA(methodName, weightFile);
        }

        void Evaluate(Particle vPositron)
        {

                SFPCAL = (vPositron.Energy(ECAL, PCAL)) / vPositron.Vector.P();
                SFECIN = (vPositron.Energy(ECAL, ECIN)) / vPositron.Vector.P();
                SFECOUT = (vPositron.Energy(ECAL, ECOUT)) / vPositron.Vector.P();

                double M2PCAL = -1.;
                double M2ECIN = -1.;
                double M2ECOUT = -1.;

                for (int i = 0; i < vPositron.Calorimeter.size(); i++)
                {
                        if (PCAL == vPositron.Calorimeter[i].layer)
                        {
                                M2PCAL = (vPositron.Calorimeter[i].m2u + vPositron.Calorimeter[i].m2v + vPositron.Calorimeter[i].m2w) / 3.;
                        }
                        if (ECIN == vPositron.Calorimeter[i].layer)
                        {
                                M2ECIN = (vPositron.Calorimeter[i].m2u + vPositron.Calorimeter[i].m2v + vPositron.Calorimeter[i].m2w) / 3.;
                        }
                        if (ECOUT == vPositron.Calorimeter[i].layer)
                        {
                                M2ECOUT = (vPositron.Calorimeter[i].m2u + vPositron.Calorimeter[i].m2v + vPositron.Calorimeter[i].m2w) / 3.;
                        }
                }

                m2PCAL = M2PCAL;
                m2ECIN = M2ECIN;
                m2ECOUT = M2ECOUT;
                
                score = PositronIdentificationReader->EvaluateMVA(methodName);
        }

        bool Accept(Particle vPositron)
        {

                //cout<<vPositron.Vector.P()<<" "<<momentumThreshold<<" "<<score<<" "<<scoreCut<<endl;
                if ((float)vPositron.Vector.P() > (float)momentumThreshold)
                        return ((float)score > (float)scoreCut);
                else return true;
        }
};

#endif
