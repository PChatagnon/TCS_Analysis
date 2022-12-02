#ifndef TCSMCEvent
#define TCSMCEvent

class MCEvent
{

public:
        // Particle and Lorentz vectors
        TLorentzVector Electron_1;
        TLorentzVector Electron_2;
        TLorentzVector Positron;
        TLorentzVector Proton;
        TLorentzVector vBeam;
        TLorentzVector vMissing;
        TLorentzVector vRestProton;
        TLorentzVector vPhoton;

        // Kinematic variables
        float t_Gen;
        float t_min_Gen;
        float MMassBeam_Gen;
        float Epho_Gen;
        float qp2_Gen;
        float M_Gen_1;
        float M_Gen_2;
        float Pt_Frac_Gen;
        float Q2_Gen;
        float real_flux_Gen;
        float virtual_flux_Gen;

        //MC weight

        float w;

        // Angular variables
        ThetaPhi cm;
        float theta_Gen;
        float phi_Gen;

        MCEvent()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_MC_Particles(hipo::bank MCEVENT, hipo::bank MCPART, bool IsGrape, bool IsJPsi)
        {

                if (IsGrape)
                {
                        Electron_1.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), me);
                        Electron_2.SetXYZM(MCPART.getFloat("px", 3), MCPART.getFloat("py", 3), MCPART.getFloat("pz", 3), me);
                        Positron.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), me);
                        Proton.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), mp);
                }

                if (!IsGrape && !IsJPsi)
                {
                        Electron_2.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), me);
                        Positron.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), me);
                        Proton.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), mp);
                        

                        float Egamma_gen = (Electron_2 + Positron + Proton - vRestProton).E();
                        float flux = n_real(ebeam, Egamma_gen);

                        float MCpsf = MCEVENT.getFloat("pbeam", 0);
			float MCcs = MCEVENT.getFloat("weight", 0);
                        // w = MCEVENT.getFloat("ebeam", 0) * MCEVENT.getFloat("weight", 0) * MCEVENT.getFloat("pbeam", 0);
                        w = MCpsf * MCcs * flux;
                }

                if (IsJPsi)
                {
                        Electron_2.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), me);
                        Positron.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), me);
                        Proton.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), mp);

                        float MC_factor_1 = MCEVENT.getFloat("ptarget", 0);
                        float MC_factor_2 = MCEVENT.getFloat("pbeam", 0);
                        float MC_factor_3 = MCEVENT.getFloat("ebeam", 0);
                        w = MC_factor_1 * MC_factor_2 * MC_factor_3;
                }

                /*if (Weighted_simu)
                {
                        w = MCEVENT.getFloat("weight", 0);
                }*/
        }

        void Get_Kinematics()
        {

                // Kinematic variables
                t_Gen = (vRestProton - Proton).M2();
                vPhoton = Proton + Positron + Electron_2 - vRestProton;
                Epho_Gen = vPhoton.E();
                vMissing = vBeam - vPhoton;
                MMassBeam_Gen = (vMissing).M2();
                qp2_Gen = (Positron + Electron_2).M2();
                M_Gen_2 = sqrt(qp2_Gen);
                M_Gen_1 = sqrt((Positron + Electron_1).M2());
                Pt_Frac_Gen = vMissing.Pt() / vMissing.P();
                Q2_Gen = 2 * ebeam * vMissing.E() * (1. - cos(vMissing.Theta()));
                t_min_Gen = T_min( 0.0, mp*mp, M_Gen_2*M_Gen_2, mp*mp, (vPhoton+vRestProton).M2());

                // Angular variables
                cm = CM(Electron_2, Positron, Proton);
                theta_Gen = cm.theta;
                phi_Gen = cm.phi;

                real_flux_Gen = n_real(ebeam, vPhoton.E());
                virtual_flux_Gen = n_virtual(ebeam, vPhoton.E(), 0.02);
        }
};

#endif
