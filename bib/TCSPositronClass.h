#ifndef TCSPositronClass
#define TCSPositronClass

class PositronIdentification
{
public:
        TString methodName;
        TString weightFile;
        double scoreCut;
        double momentumThreshold;
        Float_t P, Theta, Phi;
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

        void InitializeBDT_new_PositronIdentification()
        {
                PositronIdentificationReader->AddVariable("P", &P);
                PositronIdentificationReader->AddVariable("Theta", &Theta);
                PositronIdentificationReader->AddVariable("Phi", &Phi);
                PositronIdentificationReader->AddVariable("SFPCAL", &SFPCAL);
                PositronIdentificationReader->AddVariable("SFECIN", &SFECIN);
                PositronIdentificationReader->AddVariable("SFECOUT", &SFECOUT);
                PositronIdentificationReader->AddVariable("m2PCAL", &m2PCAL);
                PositronIdentificationReader->AddVariable("m2ECIN", &m2ECIN);
                PositronIdentificationReader->AddVariable("m2ECOUT", &m2ECOUT);
                PositronIdentificationReader->BookMVA(methodName, weightFile);
        }

        void Evaluate(Particle vLepton)
        {

                SFPCAL = (vLepton.Energy(ECAL, PCAL)) / vLepton.Vector.P();
                SFECIN = (vLepton.Energy(ECAL, ECIN)) / vLepton.Vector.P();
                SFECOUT = (vLepton.Energy(ECAL, ECOUT)) / vLepton.Vector.P();

                m2PCAL = vLepton.M2_ECAL(PCAL);
                m2ECIN = vLepton.M2_ECAL(ECIN);
                m2ECOUT = vLepton.M2_ECAL(ECOUT);
                
                score = PositronIdentificationReader->EvaluateMVA(methodName);
        }

        void Evaluate_BDT_new(Particle vLepton)
        {

                SFPCAL = (vLepton.Energy(ECAL, PCAL)) / vLepton.Vector.P();
                SFECIN = (vLepton.Energy(ECAL, ECIN)) / vLepton.Vector.P();
                SFECOUT = (vLepton.Energy(ECAL, ECOUT)) / vLepton.Vector.P();

                m2PCAL = vLepton.M2_ECAL(PCAL);
                m2ECIN = vLepton.M2_ECAL(ECIN);
                m2ECOUT = vLepton.M2_ECAL(ECOUT);

                P = vLepton.Vector.P();
                Theta = vLepton.Vector.Theta();
                Phi = vLepton.Vector.Phi();
                
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
