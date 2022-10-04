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