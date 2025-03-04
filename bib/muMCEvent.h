#ifndef ClassmuMCEvent
#define ClassmuMCEvent

class muMCEvent
{

public:
        // Particle and Lorentz vectors
        TLorentzVector Electron;
        TLorentzVector mu_plus;
        TLorentzVector mu_minus;
        TLorentzVector Proton;

        TLorentzVector vMissing;
        TLorentzVector vRestProton;
        TLorentzVector vPhoton;
        TLorentzVector vBeam;


        // Kinematic variables
        float M_Gen;
        float Q2_Gen;
        
        //MC vertex
        float vz_elec_Gen;
        float vz_mu_plus_Gen;
        float vz_mu_minus_Gen;
        float vz_prot_Gen;

        muMCEvent()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_MC_Particles(hipo::bank MCEVENT, hipo::bank MCPART, bool IsElSpectro)
        {
                if(IsElSpectro)
                {
                        mu_minus.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), mMu);
                        mu_plus.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), mMu);
                        Proton.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), mp);
                        Electron.SetXYZM(MCPART.getFloat("px", 3), MCPART.getFloat("py", 3), MCPART.getFloat("pz", 3), me);
                        
                        vz_mu_plus_Gen = MCPART.getFloat("vz", 0);
                        vz_mu_minus_Gen = MCPART.getFloat("vz", 1);
                        vz_prot_Gen = MCPART.getFloat("vz", 2);
                        vz_elec_Gen = MCPART.getFloat("vz", 3);
                }

                if(!IsElSpectro)
                {
                        mu_plus.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), mMu);
                        mu_minus.SetXYZM(MCPART.getFloat("px", 3), MCPART.getFloat("py", 3), MCPART.getFloat("pz", 3), mMu);
                        Proton.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), mp);
                        Electron.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), me);
                        
                        vz_mu_plus_Gen = MCPART.getFloat("vz", 3);
                        vz_mu_minus_Gen = MCPART.getFloat("vz", 2);
                        vz_prot_Gen = MCPART.getFloat("vz", 0);
                        vz_elec_Gen = MCPART.getFloat("vz", 1);
                }
        }

        void Get_Kinematics()
        {
                // Kinematic variables

                float qp2_Gen = (mu_plus + mu_minus).M2();
                M_Gen = sqrt(qp2_Gen);
                Q2_Gen = -1.0 * (vBeam - Electron).M2();
        }
};

#endif
