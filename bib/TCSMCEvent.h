#ifndef TCSMCEvent
#define TCSMCEvent


void GJ_angles(double &theta_GJ, double &phi_GJ, TLorentzVector vElectron, TLorentzVector vPositron, TLorentzVector vProton)
{
        TLorentzVector vBeam;
        vBeam.SetPxPyPzE(0., 0., 5, 5);

        TLorentzVector vPair;
        vPair = vPositron + vElectron;

        vPositron.Boost(-vPair.BoostVector());
        vElectron.Boost(-vPair.BoostVector());
        vProton.Boost(-vPair.BoostVector());
        vBeam.Boost(-vPair.BoostVector());

        TVector3 zAxis(0, 0, 1);
        TVector3 currentDir = vBeam.Vect().Unit(); // direction after boost
        TVector3 rotationAxis = currentDir.Cross(zAxis);
        double rotationAngle = currentDir.Angle(zAxis);

        vPositron.Rotate(rotationAngle, rotationAxis);
        vElectron.Rotate(rotationAngle, rotationAxis);
        vProton.Rotate(rotationAngle, rotationAxis);
        vBeam.Rotate(rotationAngle, rotationAxis);

        double rotation_proton = vProton.Phi();

        vPositron.Rotate(-rotation_proton, zAxis);
        vElectron.Rotate(-rotation_proton, zAxis);
        vProton.Rotate(-rotation_proton, zAxis);
        vBeam.Rotate(-rotation_proton, zAxis);

        cout << " "<< endl;
        cout << vProton.Px() << " " << vProton.Py() << " " << vProton.Pz() << " " << vProton.E() << endl;
        cout << vBeam.Px() << " " << vBeam.Py() << " " << vBeam.Pz() << " " << vBeam.E() << endl;

        theta_GJ= vElectron.Theta();
        phi_GJ=vElectron.Phi();
}



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
        float virtual_flux_Frixione_Gen;

        //MC vertex
        float vz_elec_Gen;
        float vz_posi_Gen;
        float vz_prot_Gen;

        //MC weight

        float w;

        // Angular variables
        ThetaPhi cm;
        float theta_Gen;
        float phi_Gen;
        double theta_GJ_Gen;
        double phi_GJ_Gen;

        MCEvent()
        {
                vRestProton.SetPxPyPzE(0., 0., 0., mp);
                vBeam.SetPxPyPzE(0., 0., ebeam, ebeam);
        }

        void Set_MC_Particles(hipo::bank MCEVENT, hipo::bank MCPART, bool IsGrape, bool IsJPsi, bool IsElSpectro, bool IsTCSGen)
        {

                if (IsGrape)
                {
                        Electron_1.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), me);
                        Electron_2.SetXYZM(MCPART.getFloat("px", 3), MCPART.getFloat("py", 3), MCPART.getFloat("pz", 3), me);
                        Positron.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), me);
                        Proton.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), mp);

                        vz_elec_Gen = MCPART.getFloat("vz", 3);
                        vz_posi_Gen = MCPART.getFloat("vz", 2);
                        vz_prot_Gen = MCPART.getFloat("vz", 0);
                }

                if (IsTCSGen)
                {
                        Electron_2.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), me);
                        Positron.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), me);
                        Proton.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), mp);
                        
                        vz_elec_Gen = MCPART.getFloat("vz", 0);
                        vz_posi_Gen = MCPART.getFloat("vz", 1);
                        vz_prot_Gen = MCPART.getFloat("vz", 2);

                        float Egamma_gen = (Electron_2 + Positron + Proton - vRestProton).E();
                        float flux = n_real(ebeam, Egamma_gen);

                        float MCpsf = MCEVENT.getFloat("pbeam", 0);
			float MCcs = MCEVENT.getFloat("weight", 0);
                        float MCflux = MCEVENT.getFloat("ebeam", 0);

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

                if (IsElSpectro)
                {
                        Electron_2.SetXYZM(MCPART.getFloat("px", 0), MCPART.getFloat("py", 0), MCPART.getFloat("pz", 0), me);
                        Positron.SetXYZM(MCPART.getFloat("px", 1), MCPART.getFloat("py", 1), MCPART.getFloat("pz", 1), me);
                        Proton.SetXYZM(MCPART.getFloat("px", 2), MCPART.getFloat("py", 2), MCPART.getFloat("pz", 2), mp);

                        vz_elec_Gen = MCPART.getFloat("vz", 0);
                        vz_posi_Gen = MCPART.getFloat("vz", 1);
                        vz_prot_Gen = MCPART.getFloat("vz", 2);
                        //w = 1;
                }
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

                //GJ angles
                GJ_angles(theta_GJ_Gen, phi_GJ_Gen, Electron_2,  Positron,  Proton);

                real_flux_Gen = n_real(ebeam, vPhoton.E());
                virtual_flux_Gen = n_virtual(ebeam, vPhoton.E(), 0.02);
                virtual_flux_Frixione_Gen = Frixione_wThreshold(vPhoton.E());
        }
};

#endif
