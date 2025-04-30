#ifndef ClassmuEvent
#define ClassmuEvent

class muEvent
{

public:
        // Particle and Lorentz vectors
        Particle Electron;
        Particle mu_plus;
        Particle mu_minus;
        Particle *Photons;
        TLorentzVector vBeam;
        TLorentzVector vMissing;
        TLorentzVector vRestProton;


        // Number of reconstructed particle of each type
        int np = 3; // total number of particles
        int rec_e = 0;
        int rec_mu_p = 0;
        int rec_mu_m = 0;

        // Kinematic variables
        float MMass;
        float Epho;
        float M;
        float qp2;
        float Q2;
        float t;
        float W;
        float phi;
        float theta;

        //muon ID variables
        float n_strip_PCAL_mu_plus;
	float n_strip_PCAL_mu_minus;
	float n_strip_ECIN_mu_plus;
	float n_strip_ECIN_mu_minus;
	float n_strip_ECOUT_mu_plus;	
	float n_strip_ECOUT_mu_minus;

        // run number and trigger bit
        float run;
        int trigger_bit;

        muEvent(int nb_part)
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
                np = nb_part;
                Photons = new Particle[np];
        }

        muEvent()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_nb_part(int input_np)
        {
                np = input_np;
                Photons = new Particle[np];
        }

        void Set_Particles(hipo::bank PART)
        {

                for (int i = 0; i < np; i++)
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

                        if (pid == -211)
                        {
                                mu_minus.Vector.SetXYZM(px, py, pz, mMu);
                                mu_minus.index = i;
                                mu_minus.pid = -13;
                                mu_minus.beta = beta;
                                mu_minus.status = status;
                                mu_minus.chi2 = chi2;
                                mu_minus.vertex.x = vx;
                                mu_minus.vertex.y = vy;
                                mu_minus.vertex.z = vz;
                                mu_minus.vt = vt;
                                rec_mu_m++;
                        }

                        if (pid == 211)
                        {
                                mu_plus.Vector.SetXYZM(px, py, pz, mMu);
                                mu_plus.index = i;
                                mu_plus.pid = 13;
                                mu_plus.beta = beta;
                                mu_plus.status = status;
                                mu_plus.chi2 = chi2;
                                mu_plus.vertex.x = vx;
                                mu_plus.vertex.y = vy;
                                mu_plus.vertex.z = vz;
                                mu_plus.vt = vt;
                                rec_mu_p++;
                        
                        }
                }
        }

        bool pass_topology_cut()
        {
                return (rec_mu_m == 1 && rec_mu_p == 1);
        }

        void Associate_detector_resp(hipo::bank CHE, hipo::bank SCIN, hipo::bank CALO)
        {
                vector<Particle> Particles = {Electron, mu_minus, mu_plus};

                CalorimeterResp Calo;
                CheResp Che;
                ScinResp Scin;
                for (int i = 0; i < 3; i++)
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
                                float Calohx = CALO.getFloat("hx", c);
                                float Calohy = CALO.getFloat("hy", c);
                                float Calohz = CALO.getFloat("hz", c);
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
                                        Calo.hx = Calohx;
                                        Calo.hy = Calohy;
                                        Calo.hz = Calohz;
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
                }

                Electron = Particles[0];
                mu_minus = Particles[1];
                mu_plus = Particles[2];
        }

        void Correct_Momentum(muMCEvent MC_ev)
        {
                // Correct the momentum of the muons, provided by Rafo
                TF1 *f_Eloss_mup = new TF1("f_Eloss_mup", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 80.);
                TF1 *f_Eloss_mum = new TF1("f_Eloss_mup", "[1] + [2]*(x-[0]) + [3]/(x-[0]) + [4]/((x-[0])*(x-[0]))", 0., 80.);
                f_Eloss_mup->SetParameters(-1.4651, -0.0216782, 0.00121857, -1.2927, -0.0990961);
                f_Eloss_mum->SetParameters(-1.45262, -0.0176347, 0.000494985, -1.27237, -0.111361);
                f_Eloss_mup->SetNpx(4500);
                f_Eloss_mum->SetNpx(4500);

                float corr_mup = mu_plus.Vector.P()/(f_Eloss_mup->Eval(mu_plus.Vector.P()) + 1.);
                float corr_mum = mu_minus.Vector.P()/(f_Eloss_mum->Eval(mu_minus.Vector.P()) + 1.);

                //float theta_mup = mu_plus.Vector.Theta();
                //float phi_mup = mu_plus.Vector.Phi();
                //float theta_mum = mu_minus.Vector.Theta();
                //float phi_mum = mu_minus.Vector.Phi();

                float theta_mup = MC_ev.mu_plus.Theta();
                float phi_mup = MC_ev.mu_plus.Phi();
                float theta_mum = MC_ev.mu_minus.Theta();
                float phi_mum = MC_ev.mu_minus.Phi();

                double mupPxCorr = corr_mup * sin(theta_mup) * cos(phi_mup);
		double mupPyCorr = corr_mup * sin(theta_mup) * sin(phi_mup);
		double mupPzCorr = corr_mup * cos(theta_mup);

                double mumPxCorr = corr_mum * sin(theta_mum) * cos(phi_mum);
		double mumPyCorr = corr_mum * sin(theta_mum) * sin(phi_mum);
		double mumPzCorr = corr_mum * cos(theta_mum);

                mu_plus.Vector.SetXYZM(mupPxCorr, mupPyCorr, mupPzCorr, mMu);
                mu_minus.Vector.SetXYZM(mumPxCorr, mumPyCorr, mumPzCorr, mMu);
        }


        void Get_Kinematics()
        {
                // Kinematic variables
                
                vMissing = mu_minus.Vector + mu_plus.Vector + Electron.Vector - vRestProton - vBeam;

                MMass = (vMissing).M2();
                qp2 = (mu_minus.Vector + mu_plus.Vector).M2();
                M = sqrt(qp2);
                Q2 = -1.0 * (vBeam - Electron.Vector).M2();
                t = -1.0 * (vRestProton+vMissing).M2();
                W =  (Electron.Vector - vRestProton - vBeam).M();
        }


        void Set_Run_Number(int input_run)
        {
                run = (float)(input_run);
        }


};

#endif
