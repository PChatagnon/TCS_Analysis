#ifndef TCSMomentumCorrection
#define TCSMomentumCorrection

class MomentumCorrection
{
public:
	TF1 *f1Pproton;
	TF1 *f2Pproton;
	TF1 *f3Pproton;

	TF1 *f1P;
	TF1 *f2P;
	TF1 *f3P;

	MomentumCorrection()
	{
		f1P = new TF1("f1P", "[0]+[1]*x+[2]*x*x", 0, 1);
		f1P->SetParameters(0.153319, -0.298968, 0.1607);
		f2P = new TF1("f2P", "[0]+[1]*x+[2]*x*x", 0, 1);
		f2P->SetParameters(0.0398946, -0.0748125, 0.0395764);
		f3P = new TF1("f3P", "[0]+[1]*x", 0, 1);
		f3P->SetParameters(0.0292947, -0.0577956);

		f1Pproton = new TF1("f1Pproton", "[0]+[1]*x", -180., 180.);
		f1Pproton->SetParameters(-0.00146068, -2.13735e-05);
		f2Pproton = new TF1("f2Pproton", "[0]+[1]*x", -180., 180.);
		f2Pproton->SetParameters(-0.0608671, 0.000849025);
		f3Pproton = new TF1("f3Pproton", "[0]+[1]*x", -180., 180.);
		f3Pproton->SetParameters(-0.0670748, 0.000419003);
	}

	Particle Apply_Central_Correction(Particle vProton, hipo::bank TRACK, hipo::bank TRAJ)
	{

		if (vProton.status < 4000)
			return vProton;

		Track Region1Track;
		for (int i = 0; i < TRACK.getRows(); i++)
		{

			int sector = 0;
			float chi2 = -1000.0;
			int NDF = 1;
			int indexTrack = -5000;

			if (TRACK.getInt("status", i) > 0 && TRACK.getInt("detector", i) == 5)
			{
				sector = TRACK.getInt("sector", i);
				chi2 = TRACK.getFloat("chi2", i);
				NDF = TRACK.getInt("NDF", i);
				indexTrack = TRACK.getInt("index", i);
			}

			float Nchi2 = chi2 / NDF;

			for (int c = 0; c < TRAJ.getRows(); c++)
			{

				int pindex = TRAJ.getInt("pindex", c);
				int index = TRAJ.getInt("index", c);
				int layer = TRAJ.getInt("layer", c);
				int detector = TRAJ.getInt("detector", c);
				float x = TRAJ.getFloat("x", c);
				float y = TRAJ.getFloat("y", c);
				float z = TRAJ.getFloat("z", c);

				// Region1
				if (index == (indexTrack) && detector == 5 && layer == 12 && pindex == vProton.index)
				{
					Region1Track.detector = detector;
					Region1Track.pindex = index;
					Region1Track.sector = sector;
					Region1Track.layer = layer;
					Region1Track.x = x;
					Region1Track.y = y;
					Region1Track.z = z;
					Region1Track.chi2 = Nchi2;
				}
			}
		}

		double x = Region1Track.x;
		double y = Region1Track.y;
		double z = Region1Track.z;
		double Nchi2 = Region1Track.chi2;
		float d = sqrt(x * x + y * y);
		float phi = (TMath::ATan2(y, x) * TMath::RadToDeg());
		float theta = TMath::ATan(d / z) * TMath::RadToDeg();

		double PP = vProton.Vector.P();
		double newPP = 0.0;
		if (phi > 150. || phi < -90.)
		{
			double phi1 = phi;
			if (phi > 0.)
			{
				phi1 = phi - 270.;
			}
			else
			{
				phi1 = phi + 90.;
			}
			newPP = PP * (1. - (f1Pproton->Eval(phi1)));
		}
		if (phi > -90. && phi < 30.)
		{
			newPP = PP * (1. - (f2Pproton->Eval(phi)));
		}
		if (phi < 150. && phi > 30.)
		{
			newPP = PP * (1. - (f3Pproton->Eval(phi)));
		}

		vProton.Vector.SetRho(newPP);

		return vProton;
	}

	// MC Correction

	Particle Apply_MC_Correction(Particle vProton)
	{
		double PP = vProton.Vector.P();
		double newPP = 0.0;
		if (vProton.status < 4000 && vProton.Vector.Theta() * TMath::RadToDeg() > 27.)
			newPP = PP + (f1P->Eval(PP));
		if (vProton.status < 4000 && vProton.Vector.Theta() * TMath::RadToDeg() < 27.)
			newPP = PP + (f2P->Eval(PP));
		if (vProton.status > 4000)
			newPP = PP + (f3P->Eval(PP));

		vProton.Vector.SetRho(newPP);
		return vProton;
	}
};

class Energy_loss
{
public:
	TF1 *fElectron;
	TF1 *fPositron;
	TF1 *fProtonCD;
	TF1 *fProtonFD_HighTheta;
	TF1 *fProtonFD_LowTheta;

	Energy_loss(bool inbending, bool RGA_Fall2018)
	{

		if (inbending && RGA_Fall2018)
		{

			fElectron = new TF1("fElectron", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fElectron->SetParameters(9.12163, 3.82015, -0.494202, -5.99238e-05, 2.48287e-05, -4.64859e-06);
			fPositron = new TF1("fPositron", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fPositron->SetParameters(7.41971, 3.48004, -0.571027, -4.44993e-05, 1.29072e-05, -5.74303e-06);
			fProtonCD = new TF1("fProtonCD", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonCD->SetParameters(29.1322, 4.6503, -4.38504, -9.75476e-05, 0.000154144, -3.08242e-08);
			fProtonFD_HighTheta = new TF1("fProtonFD_HighTheta", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonFD_HighTheta->SetParameters(53.038, 5.36114, -3.85085, -0.000139053, 0.000272721, -0.000193624);
			fProtonFD_LowTheta = new TF1("fProtonFD_LowTheta", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonFD_LowTheta->SetParameters(44.2175, 4.33855, -3.35995, -0.000146557, 0.000265887, -0.000174179);
		}
		else if (!inbending && RGA_Fall2018)
		{

			fElectron = new TF1("fElectron", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fElectron->SetParameters(9.12163, 3.82015, -0.494202, -5.99238e-05, 2.48287e-05, -4.64859e-06);
			fPositron = new TF1("fPositron", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fPositron->SetParameters(7.41971, 3.48004, -0.571027, -4.44993e-05, 1.29072e-05, -5.74303e-06);
			fProtonCD = new TF1("fProtonCD", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonCD->SetParameters(29.1322, 4.6503, -4.38504, -9.75476e-05, 0.000154144, -3.08242e-08);
			fProtonFD_HighTheta = new TF1("fProtonFD_HighTheta", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonFD_HighTheta->SetParameters(53.038, 5.36114, -3.85085, -0.000139053, 0.000272721, -0.000193624);
			fProtonFD_LowTheta = new TF1("fProtonFD_LowTheta", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonFD_LowTheta->SetParameters(44.2175, 4.33855, -3.35995, -0.000146557, 0.000265887, -0.000174179);
		}
		else if (!RGA_Fall2018)
		{

			fElectron = new TF1("fElectron", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fElectron->SetParameters(9.12163, 3.82015, -0.494202, -5.99238e-05, 2.48287e-05, -4.64859e-06);
			fPositron = new TF1("fPositron", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fPositron->SetParameters(7.41971, 3.48004, -0.571027, -4.44993e-05, 1.29072e-05, -5.74303e-06);
			fProtonCD = new TF1("fProtonCD", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonCD->SetParameters(29.1322, 4.6503, -4.38504, -9.75476e-05, 0.000154144, -3.08242e-08);
			fProtonFD_HighTheta = new TF1("fProtonFD_HighTheta", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonFD_HighTheta->SetParameters(53.038, 5.36114, -3.85085, -0.000139053, 0.000272721, -0.000193624);
			fProtonFD_LowTheta = new TF1("fProtonFD_LowTheta", "([0]*exp([1]+[2]*x))   * ( [5] * x * x + [4] * x + [3])");
			fProtonFD_LowTheta->SetParameters(44.2175, 4.33855, -3.35995, -0.000146557, 0.000265887, -0.000174179);
		}
	}

	void Apply_Energy_loss(Particle *Electron, Particle *Positron, Particle *Proton, int status_proton)
	{
		double corr_factor_electron = 1. - fElectron->Eval(Electron->Vector.P());
		double corr_factor_positron = 1. - fPositron->Eval(Positron->Vector.P());
		double corr_factor_proton = 1. - (Proton->Vector.Theta() * TMath::RadToDeg() < 27. ? fProtonFD_LowTheta->Eval(Proton->Vector.P()) : fProtonFD_HighTheta->Eval(Proton->Vector.P()));
		if (status_proton > 4000)
			corr_factor_proton = 1. - fProtonCD->Eval(Proton->Vector.P());
		// Electron->Vector.SetRho(Electron->Vector.P() * corr_factor_electron);
		// Positron->Vector.SetRho(Positron->Vector.P() * corr_factor_positron);
		// Proton->Vector.SetRho(Proton->Vector.P() * corr_factor_proton);

		// Electron->Vector.SetXYZM(Electron->Vector.Px()/corr_factor_electron, Electron->Vector.Py()/corr_factor_electron , Electron->Vector.Pz()/corr_factor_electron, me);
		// Positron->Vector.SetXYZM(Positron->Vector.Px()/corr_factor_positron, Positron->Vector.Py()/corr_factor_positron, Positron->Vector.Pz()/corr_factor_positron, me);
		double pPxCorr = corr_factor_proton * sin(Proton->Vector.Theta()) * cos(Proton->Vector.Phi());
		double pPyCorr = corr_factor_proton * sin(Proton->Vector.Theta()) * sin(Proton->Vector.Phi());
		double pPzCorr = corr_factor_proton * cos(Proton->Vector.Theta());
		Proton->Vector.SetXYZM(pPxCorr, pPyCorr, pPzCorr, mp);
	}
};

class Momentum_Corrections_Richard
{
public:
	std::array<double, 5> corIn_PDep_noPhiDep;
    std::array<double, 5> corOut_PDep_noPhiDep;
	Momentum_Corrections_Richard(bool inbending)
	{
		if (inbending)
		{
			corIn_PDep_noPhiDep = {9.376e-03, -7.808e-04, 1.063e-04, -2.029e-02, 4.121e-02};
			corOut_PDep_noPhiDep = {-6.520e-02, 7.099e-03, -5.929e-05, 2.145e-01, -1.153e-01};
		}

		else
		{
			corIn_PDep_noPhiDep = {-6.520e-02, 7.099e-03, -5.929e-05, 2.145e-01, -1.153e-01};
			corOut_PDep_noPhiDep = {9.376e-03, -7.808e-04, 1.063e-04, -2.029e-02, 4.121e-02};
		}
	}

	void Apply_Momentum_Corrections(Particle *Electron, Particle *Positron)
	{

		double elPCfc = Electron->Vector.P();
		double poPCfc = Positron->Vector.P();

		if (elPCfc < 1.0)
		{
			elPCfc = 1.0;
		}

		if (poPCfc < 1.0)
		{
			poPCfc = 1.0;
		}

		double elPCorr = elPCfc + (corIn_PDep_noPhiDep[0] + corIn_PDep_noPhiDep[1] * elPCfc + corIn_PDep_noPhiDep[2] * elPCfc * elPCfc + corIn_PDep_noPhiDep[3] / elPCfc + corIn_PDep_noPhiDep[4] / (elPCfc * elPCfc)) * elPCfc;
		double poPCorr = poPCfc + (corOut_PDep_noPhiDep[0] + corOut_PDep_noPhiDep[1] * poPCfc + corOut_PDep_noPhiDep[2] * poPCfc * poPCfc + corOut_PDep_noPhiDep[3] / poPCfc + corOut_PDep_noPhiDep[4] / (poPCfc * poPCfc)) * poPCfc;

		double elPxCorr = elPCorr * sin(Electron->Vector.Theta()) * cos(Electron->Vector.Phi());
		double elPyCorr = elPCorr * sin(Electron->Vector.Theta()) * sin(Electron->Vector.Phi());
		double elPzCorr = elPCorr * cos(Electron->Vector.Theta());

		Electron->Vector.SetXYZM(elPxCorr, elPyCorr, elPzCorr, me);

		double poPxCorr = poPCorr * sin(Positron->Vector.Theta()) * cos(Positron->Vector.Phi());
		double poPyCorr = poPCorr * sin(Positron->Vector.Theta()) * sin(Positron->Vector.Phi());
		double poPzCorr = poPCorr * cos(Positron->Vector.Theta());

		Positron->Vector.SetXYZM(poPxCorr, poPyCorr, poPzCorr, me);
	}

};

Particle RadiativeCorr(Particle vParticle, Particle Photons[], double thetaWin, double thetaWin1, int np)
{

	for (int i = 0; i < np; i++)
	{
		if (Photons[i].Vector.Angle(vParticle.Vector.Vect()) * TMath::RadToDeg() < thetaWin && abs(Photons[i].Vector.Theta() - vParticle.Vector.Theta()) * TMath::RadToDeg() < thetaWin1 && Photons[i].Vector.P() > 0.)
		{
			vParticle.Vector = (vParticle.Vector + Photons[i].Vector);
		}
	}
	return vParticle;
}

#endif