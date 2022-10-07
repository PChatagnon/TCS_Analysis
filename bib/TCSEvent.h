#ifndef TCSEvent
#define TCSEvent

class Event
{

public:
        // Particle and Lorentz vectors
        Particle Electron;
        Particle Positron;
        Particle Proton;
        Particle *Photons;
        TLorentzVector vBeam;
        TLorentzVector vMissing;
        TLorentzVector vRestProton;
        TLorentzVector vPhoton;

        // Event polarization
        int polarization;
        float polaT;

        // Event weight
        float w = 1;

        // Number of reconstructed particle of each type
        int np = 3; // total number of particles
        int recem = 0;
        int recep = 0;
        int recp = 0;

        // Kinematic variables
        float t;
        float MMassBeam;
        float Epho;
        float qp2;
        float M;
        float xi;
        float Pt_Frac;
        float Q2;

        // Angular variables
        ThetaPhi cm;
        float theta;
        float phi;

        // Auxiliary kinematic variables
        float b;
        float s;
        float L;
        float L0;

        // Sampling fractions
        float positron_SF;
        float electron_SF;

        // Acceptance
        float acc = 0;
        float acc_error = -1;

        // Photon flux
        float real_flux;
        float virtual_flux;

        // run number
        float run;

        Event(int nb_part)
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
                np = nb_part;
                Photons = new Particle[np];
        }

        Event()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_nb_part(int input_np)
        {
                np = input_np;
                Photons = new Particle[np];
        }

        void Set_Polarization(float pola)
        {
                polarization = pola;
        }

        void Set_Weight(float w_in)
        {
                w = w_in;
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
                        if (pid == -11)
                        {
                                Positron.Vector.SetXYZM(px, py, pz, me);
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
                                if (status > 2000)
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

                        if (pid == 22)
                        {
                                Photons[i].Vector.SetXYZM(px, py, pz, 0.);
                                Photons[i].index = i;
                                Photons[i].pid = 22;
                                Photons[i].status = status;
                                Photons[i].beta = beta;
                        }
                }
        }

        bool pass_topology_cut()
        {
                return (recem == 1 && recep == 1 && recp == 1);
        }

        void Apply_EC_Cuts(hipo::bank CALO)
        {
                Electron = ApplyECcuts(Electron, CALO);
                Positron = ApplyECcuts(Positron, CALO);
                Proton = ApplyECcuts(Proton, CALO);
        }

        bool pass_EC_cut()
        {
                return (Electron.passEC && Positron.passEC && Proton.passEC);
        }

        void Compute_SF()
        {
                electron_SF = ((Electron.Energy(ECAL, PCAL) + Electron.Energy(ECAL, ECIN) + Electron.Energy(ECAL, ECOUT))) / Electron.Vector.P();
                positron_SF = ((Positron.Energy(ECAL, PCAL) + Positron.Energy(ECAL, ECIN) + Positron.Energy(ECAL, ECOUT))) / Positron.Vector.P();
        }

        void Associate_detector_resp(hipo::bank CHE, hipo::bank SCIN)
        {
                vector<Particle> Particles = {Positron, Electron, Proton};

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
                }

                Electron = Particles[1];
                Positron = Particles[0];
                Proton = Particles[2];
        }

        void Apply_Radiative_Correction(bool is_apply_corr)
        {
                if (is_apply_corr)
                {
                        Electron = RadiativeCorr(Electron, Photons, 10., 1.5, np);
                        Positron = RadiativeCorr(Positron, Photons, 10., 1.5, np);
                }
        }

        void Apply_MC_Correction(MomentumCorrection MomCorr)
        {
                Proton = MomCorr.Apply_MC_Correction(Proton);
        }

        void Apply_Central_Correction(MomentumCorrection MomCorr, bool is_apply_corr, hipo::bank TRACK, hipo::bank TRAJ)
        {
                if (is_apply_corr)
                        Proton = MomCorr.Apply_Central_Correction(Proton, TRACK, TRAJ);
        }

        void Get_Kinematics()
        {
                // Kinematic variables
                t = (vRestProton - Proton.Vector).M2();
                vPhoton = Proton.Vector + Positron.Vector + Electron.Vector - vRestProton;
                Epho = vPhoton.E();
                vMissing = vBeam - vPhoton;
                MMassBeam = (vMissing).M2();
                qp2 = (Positron.Vector + Electron.Vector).M2();
                M = sqrt(qp2);
                Pt_Frac = vMissing.Pt() / vMissing.P();
                Q2 = 2*ebeam*vMissing.E()*(1.-cos(vMissing.Theta()));

                // Angular variables
                cm = CM(Electron, Positron, Proton);
                theta = cm.theta;
                phi = cm.phi;

                // Auxiliary kinematic variables
                double b = 2 * (Electron.Vector - Positron.Vector) * (vRestProton - Proton.Vector);
                s = (vRestProton + vPhoton).M2();
                xi = qp2 / (2 * (s - mp * mp) - qp2);
                L0 = qp2 * qp2 * sin(cm.theta * TMath::DegToRad()) * sin(cm.theta * TMath::DegToRad()) / 4.;
                L = (((qp2 - t) * (qp2 - t)) - (b * b)) / 4.;
        }

        void Get_Polarization_Transfer()
        {
                polaT = polarizationTransfer(ebeam, Epho, vPhoton.Theta());
        }

        void Set_SF(float input_electron_SF, float input_positron_SF)
        {
                electron_SF = input_electron_SF;
                positron_SF = input_positron_SF;
        }

        void Set_Kinematics(float input_t, float input_Epho, float input_qp2, float input_M,
                            float input_Pt_Frac, float input_MMassBeam, float input_theta, float input_phi, float input_s,
                            float input_xi, float input_L0, float input_L)
        {
                // Kinematic variables
                t = input_t;
                Epho = input_Epho;
                qp2 = input_qp2;
                M = input_M;
                Pt_Frac = input_Pt_Frac;
                MMassBeam = input_MMassBeam;

                // Angular variables
                theta = input_theta;
                phi = input_phi;

                // Auxiliary kinematic variables
                s = input_s;
                xi = input_xi;
                L0 = input_L0;
                L = input_L;
        }

        void Set_Vectors(TLorentzVector input_Electron, TLorentzVector input_Positron, TLorentzVector input_Proton)
        {
                Electron.Vector = input_Electron;
                Positron.Vector = input_Positron;
                Proton.Vector = input_Proton;
                vPhoton = Proton.Vector + Positron.Vector + Electron.Vector - vRestProton;
                vMissing = vBeam - vPhoton;
        }

        void Set_Acceptance(float Acc, float Acc_Error)
        {
                acc = Acc;
                acc_error = Acc_Error;
        }

        void Set_Photon_Flux(float input_real_flux, float input_virtual_flux)
        {
                real_flux = input_real_flux;
                virtual_flux = input_virtual_flux;
        }

        void Set_Run_Number(int input_run)
        {
                run = (float)(input_run);
        }
       
        void Add_Event_to_TTree()
        {
        }
};

#endif
