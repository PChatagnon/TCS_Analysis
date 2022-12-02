#ifndef TCSEEEvent
#define TCSEEEvent

class EEEvent
{

public:
        // Particle and Lorentz vectors
        Particle* Electron;
        Particle Proton;
        Particle Positron;
        TLorentzVector vBeam;
        TLorentzVector vRestProton;

        int recem = 0;
        int recp = 0;

        // Kinematic variables
        float MMassBeam;

        // Lepton HTCC Nphe
        float electron_1_Nphe = 0;
        float electron_2_Nphe = 0;

        // Lepton SF
        float electron_1_SF = 0;
        float electron_2_SF = 0;

        // run number and trigger bit
        float run;
        int trigger_bit;

        EEEvent()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
                Electron = new Particle[2];
        }

        void Set_Particles(hipo::bank PART)
        {

                for (int i = 0; i < PART.getRows(); i++)
                {

                        int pid = PART.getInt("pid", i);
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float beta = PART.getFloat("beta", i);
                        int status = abs(PART.getInt("status", i));
                        int charge = PART.getInt("charge", i);
                        float chi2 = PART.getFloat("chi2pid", i);
                        float vx = PART.getFloat("vx", i);
                        float vy = PART.getFloat("vy", i);
                        float vz = PART.getFloat("vz", i);
                        float vt = PART.getFloat("vt", i);
                        

                        if (pid == 11 && status > 2000)
                        {
                                Electron[recem].Vector.SetXYZM(px, py, pz, me);
                                Electron[recem].index = i;
                                Electron[recem].pid = 11;
                                Electron[recem].beta = beta;
                                Electron[recem].status = status;
                                Electron[recem].chi2 = chi2;
                                Electron[recem].vertex.x = vx;
                                Electron[recem].vertex.y = vy;
                                Electron[recem].vertex.z = vz;
                                Electron[recem].vt = vt;
                                recem++;
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
                }
        }

        void Compute_SF()
        {
                electron_1_SF = ((Electron[0].Energy(ECAL, PCAL) + Electron[0].Energy(ECAL, ECIN) + Electron[0].Energy(ECAL, ECOUT))) / Electron[0].Vector.P();
                electron_2_SF = ((Electron[1].Energy(ECAL, PCAL) + Electron[1].Energy(ECAL, ECIN) + Electron[1].Energy(ECAL, ECOUT))) / Electron[1].Vector.P();
        }

        void Set_Nphe_HTCC()
        {
                electron_1_Nphe = Electron[0].nphe(15);
                electron_2_Nphe = Electron[1].nphe(15);
        }

        void Apply_EC_Cuts(hipo::bank CALO)
        {
                Electron[0] = ApplyECcuts(Electron[0], CALO);
                Electron[1] = ApplyECcuts(Electron[1], CALO);
                Proton = ApplyECcuts(Proton, CALO);
        }

        void Associate_detector_resp(hipo::bank CHE, hipo::bank SCIN, hipo::bank CALO)
        {
                vector<Particle> Particles = {Positron, Electron[0], Electron[1], Proton};

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
                        Electron[0] = Particles[1];
                        Electron[1] = Particles[2];
                        Proton = Particles[3];

                }
        };
};
#endif
