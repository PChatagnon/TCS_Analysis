#ifndef TCSBGEvent
#define TCSBGEvent

class BGEvent
{

public:
        // Particle and Lorentz vectors
        Particle Electron;
        Particle Proton;
        Particle Positron;
        Particle PionM;
        TLorentzVector vBeam;
        TLorentzVector vRestProton;

        int recem = 0;
        int recep = 0;
        int recpionm = 0;
        int recp = 0;

        // Kinematic variables
        float MMassBeam;

        // Lepton HTCC Nphe
        float positron_Nphe = 0;
        float electron_Nphe = 0;
        ;

        // Lepton SF
        float positron_SF = 0;
        ;
        float electron_SF = 0;
        ;

        // run number and trigger bit
        float run;
        int trigger_bit;

        BGEvent()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        float mass_pion = 0.139;

        void Set_Particles(hipo::bank PART)
        {

                for (int i = 0; i < PART.getRows(); i++)
                {

                        int pid = PART.getInt("pid", i);
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float beta = PART.getFloat("beta", i);
                        int status = (PART.getInt("status", i));
                        int charge = PART.getInt("charge", i);
                        float chi2 = PART.getFloat("chi2pid", i);
                        float vx = PART.getFloat("vx", i);
                        float vy = PART.getFloat("vy", i);
                        float vz = PART.getFloat("vz", i);
                        float vt = PART.getFloat("vt", i);
                        if (pid == -11)
                        {
                                Positron.Vector.SetXYZM(px, py, pz, mass_pion);
                                Positron.index = i;
                                Positron.pid = -11;
                                Positron.beta = beta;
                                Positron.status = status;
                                Positron.chi2 = chi2;
                                Positron.vertex.x = vx;
                                Positron.vertex.y = vy;
                                Positron.vertex.z = vz;
                                Positron.vt = vt;
                                recep++;
                        }

                        if (pid == 11)
                        {

                                if (abs(status) > 2000)
                                {
                                        Electron.Vector.SetXYZM(px, py, pz, me);
                                        Electron.index = i;
                                        Electron.pid = 11;
                                        Electron.beta = beta;
                                        Electron.status = status;
                                        Electron.chi2 = chi2;
                                        Electron.vertex.x = vx;
                                        Electron.vertex.y = vy;
                                        Electron.vertex.z = vz;
                                        Electron.vt = vt;
                                        recem++;
                                }
                        }

                        if (pid == 2212)
                        {
                                Proton.Vector.SetXYZM(px, py, pz, mp);
                                Proton.index = i;
                                Proton.pid = 2212;
                                Proton.beta = beta;
                                Proton.status = status;
                                Proton.chi2 = chi2;
                                Proton.vertex.x = vx;
                                Proton.vertex.y = vy;
                                Proton.vertex.z = vz;
                                Proton.vt = vt;
                                recp++;
                        }

                        if (pid == -211)
                        {
                                PionM.Vector.SetXYZM(px, py, pz, mass_pion);
                                PionM.index = i;
                                PionM.pid = -211;
                                PionM.status = status;
                                PionM.beta = beta;
                                recpionm++;
                        }
                }
        }

        void Compute_SF()
        {
                if (abs(Electron.status) > 2000)
                        electron_SF = ((Electron.Energy(ECAL, PCAL) + Electron.Energy(ECAL, ECIN) + Electron.Energy(ECAL, ECOUT))) / Electron.Vector.P();
                positron_SF = ((Positron.Energy(ECAL, PCAL) + Positron.Energy(ECAL, ECIN) + Positron.Energy(ECAL, ECOUT))) / Positron.Vector.P();
        }

        void Set_Nphe_HTCC()
        {
                positron_Nphe = Positron.nphe(15);
                electron_Nphe = Electron.nphe(15);
        }

        void Apply_EC_Cuts(hipo::bank CALO)
        {
                Electron = ApplyECcuts(Electron, CALO);
                Positron = ApplyECcuts(Positron, CALO);
                Proton = ApplyECcuts(Proton, CALO);
                PionM = ApplyECcuts(PionM, CALO);
        }

        bool pass_EC_cut()
        {
                return (Electron.passEC && Positron.passEC && Proton.passEC && PionM.passEC);
        }

        void Associate_detector_resp(hipo::bank CHE, hipo::bank SCIN, hipo::bank CALO)
        {
                vector<Particle> Particles = {Positron, Electron, Proton, PionM};

                CalorimeterResp Calo;
                CheResp Che;
                ScinResp Scin;
                for (int i = 0; i < Particles.size(); i++)
                {

                        for (int c = 0; c < CHE.getRows(); c++)
                        {
                                int Chepindex = CHE.getInt("pindex", c);
                                int Chedetector = CHE.getInt("detector", c);
                                int Chesector = CHE.getInt("sector", c);
                                float Chenphe = CHE.getFloat("nphe", c);
                                float Chetime = CHE.getFloat("time", c);
                                float Chechi2 = CHE.getFloat("chi2", c);
                                float Chex = CHE.getFloat("x", c);
                                float Chey = CHE.getFloat("y", c);
                                float Chez = CHE.getFloat("z", c);

                                if (Chepindex == (Particles[i].index))
                                {
                                        Che.detector = Chedetector;
                                        Che.pindex = Chepindex;
                                        Che.sector = Chesector;
                                        Che.nphe = Chenphe;
                                        Che.time = Chetime;
                                        Che.chi2 = Chechi2;
                                        Che.x = Chex;
                                        Che.y = Chey;
                                        Che.z = Chez;
                                        Particles[i].Cherenkov.push_back(Che);
                                }
                        }
                        for (int c = 0; c < SCIN.getRows(); c++)
                        {
                                int Scindetector = SCIN.getInt("detector", c);
                                int Scinpindex = SCIN.getInt("pindex", c);
                                float Scintime = SCIN.getFloat("time", c);
                                float Scinpath = SCIN.getFloat("path", c);
                                float Scinenergy = SCIN.getFloat("energy", c);
                                int Scinsector = SCIN.getInt("sector", c);

                                if (Scinpindex == (Particles[i].index))
                                {
                                        Scin.detector = Scindetector;
                                        Scin.pindex = Scinpindex;
                                        Scin.t = Scintime;
                                        Scin.path = Scinpath;
                                        Scin.energy = Scinenergy;
                                        Scin.sector = Scinsector;

                                        if (Particles[i].Scintillator.energy < Scinenergy)
                                        {
                                                Particles[i].Scintillator = Scin;
                                        };
                                }
                        }

                        for (int c = 0; c < CALO.getRows(); c++)
                        {
                                int Calopindex = CALO.getInt("pindex", c);
                                int Calosector = CALO.getInt("sector", c);
                                int Calolayer = CALO.getInt("layer", c);
                                int Calodetector = CALO.getInt("detector", c);
                                float Caloenergy = CALO.getFloat("energy", c);
                                float Calox = CALO.getFloat("x", c);
                                float Caloy = CALO.getFloat("y", c);
                                float Caloz = CALO.getFloat("z", c);
                                float Calou = CALO.getFloat("lu", c);
                                float Calov = CALO.getFloat("lv", c);
                                float Calow = CALO.getFloat("lw", c);
                                float Calodu = CALO.getFloat("du", c);
                                float Calodv = CALO.getFloat("dv", c);
                                float Calodw = CALO.getFloat("dw", c);
                                float Calom2u = CALO.getFloat("m2u", c);
                                float Calom2v = CALO.getFloat("m2v", c);
                                float Calom2w = CALO.getFloat("m2w", c);
                                float Calom3u = CALO.getFloat("m3u", c);
                                float Calom3v = CALO.getFloat("m3v", c);
                                float Calom3w = CALO.getFloat("m3w", c);

                                if (Calopindex == (Particles[i].index))
                                {
                                        Calo.detector = Calodetector;
                                        Calo.pindex = Calopindex;
                                        Calo.sector = Calosector;
                                        Calo.layer = Calolayer;
                                        Calo.energy = Caloenergy;
                                        Calo.x = Calox;
                                        Calo.y = Caloy;
                                        Calo.z = Caloz;
                                        Calo.u = Calou;
                                        Calo.v = Calov;
                                        Calo.w = Calow;
                                        Calo.du = Calodu;
                                        Calo.dv = Calodv;
                                        Calo.dw = Calodw;
                                        Calo.m2u = Calom2u;
                                        Calo.m2v = Calom2v;
                                        Calo.m2w = Calom2w;
                                        Calo.m3u = Calom3u;
                                        Calo.m3v = Calom3v;
                                        Calo.m3w = Calom3w;

                                        Particles[i].Calorimeter.push_back(Calo);
                                }
                        }

                        Positron = Particles[0];
                        Electron = Particles[1];
                        Proton = Particles[2];
                        PionM = Particles[3];
                }
        };
};
#endif
