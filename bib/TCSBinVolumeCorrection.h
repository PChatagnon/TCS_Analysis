#ifndef TCSBinVolumeCorrection
#define TCSBinVolumeCorrection

class BinVolumeCorrection
{
public:
        double *volume1;
        double *volume2;

        double limitminPhi1 = -40.;
        double limitmaxPhi1 = 40.;
        double limitminTheta1 = 50.;
        double limitmaxTheta1 = 80.;

        double limitminPhi2 = 140.;
        double limitmaxPhi2 = 220.;
        double limitminTheta2 = 100.;
        double limitmaxTheta2 = 130.;

        BinVolumeCorrection(Acceptance input_Acc)
        {
                int bin_t=input_Acc.bin_t;
                int bin_M=input_Acc.bin_M;
                int bin_Eg=input_Acc.bin_Eg;

                volume1 = new double[bin_t * bin_M * bin_Eg];
                volume2 = new double[bin_t * bin_M * bin_Eg];

                for (int Egbin = 0; Egbin < bin_Eg; Egbin++)
                {
                        for (int mass = 0; mass < bin_M; mass++)
                        {
                                for (int tbinning = 0; tbinning < bin_t; tbinning++)
                                {

                                        double vol1 = 0.;
                                        double vol2 = 0.;
                                        for (int phi = 0; phi < 100; phi++)
                                        {
                                                for (int theta = 0; theta < 100; theta++)
                                                {
                                                        double phiValue1 = limitminPhi1 + ((limitmaxPhi1 - limitminPhi1) / 99.) * phi;
                                                        double thetaValue1 = limitminTheta1 + ((limitmaxTheta1 - limitminTheta1) / 99.) * theta;
                                                        double phiValue2 = limitminPhi2 + ((limitmaxPhi2 - limitminPhi2) / 99.) * phi;
                                                        if (phiValue2 > 180.)
                                                                phiValue2 = phiValue2 - 360.;
                                                        double thetaValue2 = limitminTheta2 + ((limitmaxTheta2 - limitminTheta2) / 99.) * theta;

                                                        double AccValue1 = input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetBinContent((input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetXaxis()->FindBin(phiValue1)), (input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetYaxis()->FindBin(thetaValue1)));
                                                        double AccValue2 = input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetBinContent((input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetXaxis()->FindBin(phiValue2)), (input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetYaxis()->FindBin(thetaValue2)));

                                                        double ErrAccValue1 = input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetBinError((input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetXaxis()->FindBin(phiValue1)), (input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetYaxis()->FindBin(thetaValue1)));
                                                        double ErrAccValue2 = input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetBinError((input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetXaxis()->FindBin(phiValue2)), (input_Acc.Acc[tbinning + bin_t * mass + bin_t * bin_M * Egbin]->GetYaxis()->FindBin(thetaValue2)));
                                                        if (AccValue1 > 0.05 && (ErrAccValue1 / AccValue1) < 0.5)
                                                                vol1 = vol1 + 1.;
                                                        if (AccValue2 > 0.05 && (ErrAccValue2 / AccValue2) < 0.5)
                                                                vol2 = vol2 + 1.;
                                                }
                                        }

                                        volume1[tbinning + bin_t * mass + bin_t * bin_M * Egbin] = (vol1 / 10000.);
                                        volume2[tbinning + bin_t * mass + bin_t * bin_M * Egbin] = (vol2 / 10000.);
                                }
                        }
                }
        }

       /* Print(){
            for (int Egbin = 0; Egbin < Acc_TCS.bin_Eg; Egbin++)
	{
		for (int mass = 0; mass < Acc_TCS.bin_M; mass++)
		{
			for (int tbinning = 0; tbinning < Acc_TCS.bin_t; tbinning++)
			{
				cout << "volume 1 " << Egbin << " " << mass << " " << tbinning << " " << volume1[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin] << " " << endl;
				cout << "volume 2 " << Egbin << " " << mass << " " << tbinning << " " << volume2[tbinning + Acc_TCS.bin_t * mass + Acc_TCS.bin_t * Acc_TCS.bin_M * Egbin] << " " << endl;
			}
		}
	}    
        }*/
};

#endif
