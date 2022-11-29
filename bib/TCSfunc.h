#ifndef TCSfunc
#define TCSfunc

ThetaPhi CM(TLorentzVector vElectron, TLorentzVector vPositron, TLorentzVector vProton)
{

	TLorentzVector vRestProton;
	TLorentzVector vPhoton;
	TLorentzVector vProtonLF;
	TLorentzVector vElectronLF;
	TLorentzVector vPositronLF;
	TLorentzVector vRestProtonLF;
	TLorentzVector vPhotonLF;
	TLorentzVector vLepton;
	TLorentzVector vIncoming;

	ThetaPhi CM;

	double phiCM;
	double thetaCM;
	double Pi = 3.14159265359;

	TLorentzVector vBeam;
	vBeam.SetPxPyPzE(0, 0, 10.6, 10.6);
	vRestProton.SetPxPyPzE(0, 0, 0, 0.938);
	vPhoton = vProton + vPositron + vElectron - vRestProton;
	phiCM = 0.0;
	thetaCM = 0.0;

	vIncoming = vProton + vPositron + vElectron;
	vPhoton.Boost(-vIncoming.BoostVector());
	vRestProton.Boost(-vIncoming.BoostVector());
	vProton.Boost(-vIncoming.BoostVector());
	vElectron.Boost(-vIncoming.BoostVector());
	vPositron.Boost(-vIncoming.BoostVector());
	if (vProton.X() == 0 && vProton.Y() == 0 && vProton.Z() == 0)
		cout << "problem1" << endl;

	vPhoton.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vElectron.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vPositron.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vRestProton.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	vProton.Rotate(Pi + vProton.Angle(vBeam.Vect().Unit()), vProton.Vect().Cross(vBeam.Vect().Unit()));
	// cout<<"problem proton1"<<endl;
	vElectron.RotateZ(-vRestProton.Phi());
	vPositron.RotateZ(-vRestProton.Phi());
	vProton.RotateZ(-vRestProton.Phi());
	vPhoton.RotateZ(-vRestProton.Phi());
	vRestProton.RotateZ(-vRestProton.Phi());

	phiCM = (vElectron.Phi());

	vLepton = vPositron + vElectron;
	vElectron.Boost(-vLepton.BoostVector());
	vProton.Boost(-vLepton.BoostVector());
	vPhoton.Boost(-vLepton.BoostVector());
	vRestProton.Boost(-vLepton.BoostVector());
	vPositron.Boost(-vLepton.BoostVector());

	thetaCM = vElectron.Angle(-vProton.Vect().Unit());

	CM.theta = thetaCM * TMath::RadToDeg();
	CM.phi = phiCM * TMath::RadToDeg();
	return CM;
}

ThetaPhi CM(Particle vElectron, Particle vPositron, Particle vProton)
{

	TLorentzVector vRestProton;
	TLorentzVector vPhoton;
	TLorentzVector vProtonLF;
	TLorentzVector vElectronLF;
	TLorentzVector vPositronLF;
	TLorentzVector vRestProtonLF;
	TLorentzVector vPhotonLF;
	TLorentzVector vLepton;
	TLorentzVector vIncoming;

	ThetaPhi CM;

	double phiCM;
	double thetaCM;
	double Pi = 3.14159265359;

	TLorentzVector vBeam;
	vBeam.SetPxPyPzE(0, 0, 10.6, 10.6);
	vRestProton.SetPxPyPzE(0, 0, 0, 0.938);
	vPhoton = vProton.Vector + vPositron.Vector + vElectron.Vector - vRestProton;
	phiCM = 0.0;
	thetaCM = 0.0;

	vIncoming = vProton.Vector + vPositron.Vector + vElectron.Vector;
	vPhoton.Boost(-vIncoming.BoostVector());
	vRestProton.Boost(-vIncoming.BoostVector());
	vProton.Vector.Boost(-vIncoming.BoostVector());
	vElectron.Vector.Boost(-vIncoming.BoostVector());
	vPositron.Vector.Boost(-vIncoming.BoostVector());
	if (vProton.Vector.X() == 0 && vProton.Vector.Y() == 0 && vProton.Vector.Z() == 0)
		cout << "problem1" << endl;

	// cout<<"problem proton"<<endl;
	// cout<<"Incoming "<< vIncoming.X()<<" "<<vIncoming.Y()<<" "<<vIncoming.Z()<<" "<<endl;
	// cout<<"proton "<< vProton.Vector.X()<<" "<<vProton.Vector.Y()<<" "<<vProton.Vector.Z()<<" "<<vProton.Vector.M()<<" "<<endl;
	vPhoton.Rotate(Pi + vProton.Vector.Angle(vBeam.Vect().Unit()), vProton.Vector.Vect().Cross(vBeam.Vect().Unit()));
	vElectron.Vector.Rotate(Pi + vProton.Vector.Angle(vBeam.Vect().Unit()), vProton.Vector.Vect().Cross(vBeam.Vect().Unit()));
	vPositron.Vector.Rotate(Pi + vProton.Vector.Angle(vBeam.Vect().Unit()), vProton.Vector.Vect().Cross(vBeam.Vect().Unit()));
	vRestProton.Rotate(Pi + vProton.Vector.Angle(vBeam.Vect().Unit()), vProton.Vector.Vect().Cross(vBeam.Vect().Unit()));
	vProton.Vector.Rotate(Pi + vProton.Vector.Angle(vBeam.Vect().Unit()), vProton.Vector.Vect().Cross(vBeam.Vect().Unit()));
	// cout<<"problem proton1"<<endl;
	vElectron.Vector.RotateZ(-vRestProton.Phi());
	vPositron.Vector.RotateZ(-vRestProton.Phi());
	vProton.Vector.RotateZ(-vRestProton.Phi());
	vPhoton.RotateZ(-vRestProton.Phi());
	vRestProton.RotateZ(-vRestProton.Phi());

	phiCM = (vElectron.Vector.Phi());

	vLepton = vPositron.Vector + vElectron.Vector;
	vElectron.Vector.Boost(-vLepton.BoostVector());
	vProton.Vector.Boost(-vLepton.BoostVector());
	vPhoton.Boost(-vLepton.BoostVector());
	vRestProton.Boost(-vLepton.BoostVector());
	vPositron.Vector.Boost(-vLepton.BoostVector());

	thetaCM = vElectron.Vector.Angle(-vProton.Vector.Vect().Unit());

	CM.theta = thetaCM * TMath::RadToDeg();
	CM.phi = phiCM * TMath::RadToDeg();
	return CM;
}
int bint(double t, double limit[3])
{

	if (-t < limit[0])
		return 0;
	if (-t < limit[1])
		return 1;
	if (-t < limit[2])
		return 2;
	if (-t > limit[2])
		return 3;
	else
		return -1;
}

int binP(double P)
{

	double Pbins[3] = {2., 3., 4.};
	if (P < Pbins[0])
		return 0;
	if (P < Pbins[1])
		return 1;
	if (P < Pbins[2])
		return 2;
	if (P > Pbins[2])
		return 3;
	else
		return -1;
}

double Sigma(vector<double> values)
{
	double ave = 0.0;
	double sigma = 0.0;
	int nb = 0;
	for (int i = 0; i < values.size(); i++)
	{
		if (values[i] != 0.0)
		{
			nb++;
			ave += values[i];
		}
	}
	ave = ave / nb;

	for (int i = 0; i < values.size(); i++)
	{
		if (values[i] != 0.0)
		{
			sigma += ((values[i] - ave) * (values[i] - ave));
		}
	}
	sigma = sigma / nb;
	sigma = sqrt(sigma);

	return sigma;
}

double VT(double mom, double mass, double time, double path, double STT)
{

	float beta = mom / sqrt(mom * mom + mass * mass);
	float vt = time - STT - path / 29.92f / beta;
	return vt;
}
bool goodElectron(Particle p)
{

	return (p.status > 2000);
}

bool goodProton(Particle p)
{

	return true;
}

TGraphErrors *MaxAcc(double cut, TH2D *Acc)
{

	THStack *hs1 = new THStack(Acc, "Y");
	TGraphErrors *AccLimitMax = new TGraphErrors(); // = TGraphErrors();
	AccLimitMax->SetMarkerColor(2);

	for (int j = 0; j < 10; j++)
	{
		TH1F *slice = (TH1F *)hs1->GetHists()->At(j);
		double maxAcc = 180.0;
		// look for maximum of fiducial acceptance
		for (int l = 20; l > 0; l--)
		{
			double value = slice->GetBinContent(l);
			// double error=slice->GetBinError(l);
			if (value > cut)
			{
				maxAcc = slice->GetBinCenter(l + 1);
				break;
			}
		}

		AccLimitMax->SetPoint(j, Acc->GetXaxis()->GetBinCenter(j + 1), maxAcc);
		// cout<<maxAcc<<endl;
	}

	return AccLimitMax;
}

TGraphErrors *MinAcc(double cut, TH2D *Acc)
{

	THStack *hs1 = new THStack(Acc, "Y");
	TGraphErrors *AccLimitMin = new TGraphErrors(); // = TGraphErrors();
	AccLimitMin->SetMarkerColor(4);

	for (int j = 0; j < 10; j++)
	{
		TH1F *slice = (TH1F *)hs1->GetHists()->At(j);
		double minAcc = 0.0;
		// look for maximum of fiducial acceptance
		for (int k = 1; k < 20 + 1; k++)
		{
			double value = slice->GetBinContent(k);
			// double error=slice->GetXaxis()->GetBinError(k);
			if (value > cut)
			{
				minAcc = slice->GetBinCenter(k - 1);
				break;
			}
		}
		AccLimitMin->SetPoint(j, Acc->GetXaxis()->GetBinCenter(j + 1), minAcc);
	}

	return AccLimitMin;
}

double smear(TH1F histoIn, int bin)
{

	double maxrange = histoIn.GetXaxis()->GetXmax();
	double minrange = histoIn.GetXaxis()->GetXmin();
	int binNum = (histoIn.GetSize());
	int nstepSmearing = 10000;
	TH1F *patate = new TH1F(Form("amp%d", bin), ";R^{R};Counts", 300, -1., 1.);
	cout << "start smear " << endl;
	for (int it = 0; it < nstepSmearing; it++)
	{
		TRandom random = TRandom(it);

		double Rc = 0.0;
		double R0 = 0.0;
		double R = 0.0;

		for (int i = 1; i < binNum - 1; i++)
		{
			double meanY = histoIn.GetBinContent(i);
			double phi = histoIn.GetBinCenter(i);
			// cout<<"phi "<<phi<<" mean "<<meanY<<endl;
			double sigmaY = histoIn.GetBinError(i);
			double valueY = random.Gaus(meanY, sigmaY);

			Rc += valueY * cos(phi * TMath::DegToRad());
			R0 += valueY;
		}

		// cout<<bin<<" "<<Rc<<" "<<R0<<endl;

		R = Rc / R0;

		patate->Fill(R);
		// cout<<R<<endl;
	}

	gStyle->SetOptStat(1111);
	TCanvas *patatec = new TCanvas("aaddr ", "A", 3000, 2000);
	// patate->Draw();
	patate->Fit("gaus"); //,"q0"
	patatec->SaveAs(Form("amp%d.pdf", bin));
	gStyle->SetOptStat(111);
	cout << "end smear " << endl;
	cout << patate->GetFunction("gaus")->GetParameter(1) << endl;
	cout << patate->GetFunction("gaus")->GetParameter(2) << endl;
	return patate->GetFunction("gaus")->GetParameter(2);
}

double smear(TH1F histoIn, TH1F EhistoIn, int bin)
{

	double maxrange = histoIn.GetXaxis()->GetXmax();
	double minrange = histoIn.GetXaxis()->GetXmin();
	int binNum = (histoIn.GetSize());
	int nstepSmearing = 10;
	TH1F *patate = new TH1F(Form("amp%d", bin), "amplitude;R^{R};Counts", 300, -2., 2.);
	cout << "start smear " << endl;
	for (int it = 0; it < nstepSmearing; it++)
	{
		TRandom random = TRandom(it);

		double Rc = 0.0;
		double R0 = 0.0;
		double R = 0.0;

		for (int i = 1; i < binNum - 1; i++)
		{
			double meanY = histoIn.GetBinContent(i);
			double phi = histoIn.GetBinCenter(i);
			double sigmaY = EhistoIn.GetBinContent(i);
			double valueY = random.Gaus(meanY, sigmaY);

			Rc += valueY * cos(phi * TMath::DegToRad());
			R0 += valueY;
		}

		cout << bin << " " << Rc << " " << R0 << endl;

		R = Rc / R0;

		patate->Fill(R);
		// cout<<R<<endl;
	}

	patate->Fit("gaus", "q0");

	TCanvas *patatec = new TCanvas("aaddr ", "A", 4000, 2000);
	patate->Draw();
	patatec->SaveAs(Form("amp%d.pdf", bin));

	cout << "end smear " << endl;
	cout << patate->GetFunction("gaus")->GetParameter(1) << endl;
	cout << patate->GetFunction("gaus")->GetParameter(2) << endl;
	return patate->GetFunction("gaus")->GetParameter(2);
}

PairValue smearFit(TH1F histoIn, int bin)
{

	int binNum = (histoIn.GetSize());
	int nstepSmearing = 10000;
	TH1F *patate = new TH1F(Form("amp fit %d", bin), "amplitude;Fit;Counts", 1000, -1.0, 0.5);
	cout << "start smear fit" << endl;

	PairValue results;
	for (int it = 0; it < nstepSmearing; it++)
	{
		TRandom random = TRandom(it);

		TH1F *fitTemp = new TH1F(Form("temp%d%d", it, bin), Form("temp%d", it), 10, -180, 180);
		for (int i = 1; i < binNum - 1; i++)
		{
			double meanY = histoIn.GetBinContent(i);
			double sigmaY = histoIn.GetBinError(i);

			double valueY = random.Gaus(meanY, sigmaY);
			fitTemp->SetBinContent(i, valueY);
		}

		TF1 *f2 = new TF1("f2", "[0]*sin(x*[1]/[2])", -180, 180);
		f2->FixParameter(1, 3.14159264);
		f2->FixParameter(2, 180);
		fitTemp->Fit("f2", "q0"); //->Draw();//

		double fitresult = fitTemp->GetFunction("f2")->GetParameter(0);
		double fiterror = fitTemp->GetFunction("f2")->GetParError(0);

		patate->Fill(fitresult);
	}

	gStyle->SetOptStat(1111);
	TCanvas *patatec = new TCanvas("aaddr ", "A", 4000, 2000);
	// patate->Draw();
	patate->Draw(); // Fit("gaus","q");
	patatec->SaveAs(Form("ampfit%d.pdf", bin));
	gStyle->SetOptStat(111);
	cout << "end smear " << endl;
	cout << patate->GetMean() << endl;
	cout << patate->GetStdDev() << endl;
	results.mean = patate->GetMean();	 //->GetFunction("gaus")->GetParameter(1);
	results.sigma = patate->GetStdDev(); //->GetFunction("gaus")->GetParameter(2);
	return results;
}

double polarizationTransfer(double Eb, double Eg, double angle)
{

	double u = angle * Eb;
	double xi = 1. / (1. + u * u);
	double E1 = Eb;
	double k = Eg;
	double E2 = Eb - Eg;
	double delta = k / (2. * E1 * E2);
	double BigDelta = (6. / 121.) * (xi / delta);

	double x_tab[19] = {0.5, 1., 2., 4., 8., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 120.};
	double y_tab[19] = {0.0145, 0.0490, 0.14, 0.3312, 0.6758, 1.126, 1.367, 1.564, 1.731, 1.875, 2.001, 2.114, 2.216, 2.393, 2.545, 2.676, 2.793, 2.897, 3.078};
	auto *F = new TGraph(19, x_tab, y_tab);

	double evalF = F->Eval(BigDelta);

	double a = (1. / 137.);
	double fZ = a * a * ((1. / (1. + a * a)) + 0.20206 - 0.0369 * a * a + 0.0083 * a * a * a * a - 0.002 * a * a * a * a * a * a);

	double gamma = evalF - log(delta) - 2. - fZ;

	double I = (E1 * E1 + E2 * E2) * (3. + 2. * gamma) - 2. * E1 * E2 * (1. + 4. * u * u * xi * xi * gamma);
	double L = k * ((E1 + E2) * (3. + 2. * gamma) - 2. * E2 * (1. + 4. * u * u * xi * xi * gamma)) / I;

	double Pola = L;
	if (I < 0.01 || Pola > 1.)
		Pola = 1.;

	return Pola;
}

double n_real(double Eb, double Eg)
{
	return 0.5 * (5.0 / 929.0) * (1 / Eg) * ((4.0 / 3.0) - (4.0 / 3.0) * (Eg / Eb) + (Eg * Eg) / (Eb * Eb));
}

double n_virtual(double Eb, double Eg, double Q2_max)
{
	double alpha = 1. / 137.;
	double PI = 3.14159265358979312;
	double x = Eg / Eb;
	double me = 0.00051;
	double Mp = 0.9383;
	double Q2_min = me * me * x * x / (1 - x);
	return (1 / Eb) * alpha / (PI * x) * ((1 - x + x * x / 2) * log(Q2_max / Q2_min) - (1 - x));
}

double Lambda( double x, double y, double z )
{
  return (x - y - z)*(x - y - z) - 4*y*z;
}

//From Byukling Kayanti Formula (5.14) Page 86
double T_min( double ma_2, double mb_2, double m1_2, double m2_2, double s) // arguments are squares of masses of particles in the reaction a+b->1+2, and s is the square of the total c.m. energy i.e. (a+b)^2
{
  return ma_2 + m1_2 - (1/(2*s))*( (s + ma_2 - mb_2)*(s + m1_2 - m2_2) - sqrt( Lambda(s, ma_2, mb_2)*Lambda(s, m1_2, m2_2) ) );
}

#endif
