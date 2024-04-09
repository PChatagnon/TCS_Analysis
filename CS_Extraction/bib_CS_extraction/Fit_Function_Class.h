#ifndef Fit_Function_class
#define Fit_Function_class

double crystalball_function(double x, double mean, double sigma, double alpha, double n) {
        // evaluate the crystal ball function
        if (sigma < 0.)     
                return 0.;
        double z = (x - mean)/sigma; 
        if (alpha < 0)
                z = -z; 
        double abs_alpha = std::abs(alpha);
        double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
        double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
        double N = 1./(sigma*(C+D));
        if (z  > - abs_alpha)
                return std::exp(- 0.5 * z * z);
        else {
                //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
                double nDivAlpha = n/abs_alpha;
                double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
                double B = nDivAlpha -abs_alpha;
                double arg = nDivAlpha/(B-z);
                return AA * std::pow(arg,n);
        }
}

class Fit_Function
{
public:
        TF1 *function;
        TF1 *function_Signal;
        TF1 *function_BG;
        float min_fit;
        float max_fit;
        TH1D *input_Data_hist;
        double chi2;
        double NDF;
        TMatrixDSym cov_matrix;
        int nb_param;
        // TFitResultPtr fitResult;

        Fit_Function() {}

        void Set_Data_hist(TH1D *in_Data_hist)
        {
                input_Data_hist = (TH1D *)in_Data_hist->Clone("Data_hist");
        }

        void Set_Limits(float in_min_fit, float in_max_fit)
        {
                min_fit = in_min_fit;
                max_fit = in_max_fit;
        }

        void Draw_Functions()
        {
                cout << "Draw function" << endl;
                //this->function->SetFillColorAlpha(kRed, 0.8);
                //this->function->Draw("3 same"); // "E3" for a one-sigma band
                this->function->Draw("same l");
                cout << "Draw BG" << endl;
                this->function_BG->Draw("same l");
                cout << "Draw signal" << endl;
                this->function_Signal->Draw("same l");
                cout << "End drawing" << endl;
        }

        float Get_Integral_Signal()
        {
                cout<<"NB jpsi"<<endl;
                cout<<function_Signal->GetParameter(0)/(double)(input_Data_hist->GetXaxis()->GetBinWidth(2))<<endl;
                cout<<(function_Signal->Integral(0., 10.)) / (input_Data_hist->GetXaxis()->GetBinWidth(2))<<endl;
                cout<<function_Signal->GetParameter(0)<<endl;
                return function_Signal->GetParameter(0)/(double)(input_Data_hist->GetXaxis()->GetBinWidth(2));


                /*return (function_Signal->Integral(0., 10.)) / (input_Data_hist->GetXaxis()->GetBinWidth(2));*/
        }

        float Get_Integral_Error_Signal()
        {
                /*cout << "IN COVARIANCE MATRICES" << endl;
                double epsilon = 1.e-7;
                TMatrixDSym covGG = cov_matrix.GetSub(0, function_Signal->GetNpar() - 1, 0, function_Signal->GetNpar() - 1); //= cov_matrix.GetSub(0, , 0, function_Signal->GetNpar());
                covGG.Print();
                double errNgg = function_Signal->IntegralError(0., 10., function_Signal->GetParameters(), covGG.GetMatrixArray(),epsilon) / (double)(input_Data_hist->GetXaxis()->GetBinWidth(2));
                cout << "IN ERROR CAL" << endl;
                cout << "integral " << function_Signal->Integral(0., 10.) << endl;
                cout << "error integral " << function_Signal->IntegralError(0., 10., function_Signal->GetParameters(), covGG.GetMatrixArray(),epsilon) << endl;
                cout << errNgg << " " << sqrt(this->Get_Integral_Signal()) << endl;
                double result = TMath::Max(errNgg, (double)sqrt(this->Get_Integral_Signal()));
                return errNgg; */

                cout<<"error NB jpsi"<<endl;
                cout<<function->GetParError(0)/(double)(input_Data_hist->GetXaxis()->GetBinWidth(2))<<endl;
                cout<<function->GetParError(0)<<endl;
                return function->GetParError(0)/(double)(input_Data_hist->GetXaxis()->GetBinWidth(2));

        }

        void Crystall_Ball_fit(TString options, TString name){

                int nb_param = 8;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                cout << "name BG func" << endl;
                cout << name_bg << endl;

                
                //function = new TF1(name_function,"[0]*crystalball_function(x, [1], [2], [3], [4]) + ([5]+[6]*(x-3.1) + [7]*(x-3.1)*(x-3.1))",min_fit, max_fit);
                function = new TF1(name_function,"[0]*ROOT::Math::crystalball_function(x,[3], [4], [2], [1]) + ([5]+[6]*(x-3.1) + [7]*(x-3.1)*(x-3.1))",min_fit, max_fit);
                
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100.);
                function->SetParameter (1,3.096);
                function->SetParLimits (1,3.02, 3.2);
                function->SetParameter (2,0.04);
                function->SetParLimits (2,0.025, 0.15);
                function->SetParameter(3,0.75);
                //function->SetParLimits(3,0.0,1.5);
                function->SetParameter(4,500);
                function->SetParLimits (4,1,5000);
                //BG
                function->SetParameter(5, 7.);
                function->SetParLimits(5, 0.00, 100000.);
                function->SetParameter(6, -2.);
                //function->SetParLimits(4, -100., 100);
                function->SetParameter(7, 4.);
                function->SetParLimits(7, 0., 100000.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();

                function_Signal = new TF1(name_signal,"[0]*crystalball_function(x, [1], [2], [3], [4])", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetParameter(3, function->GetParameter(3));
                function_Signal->SetParameter(4, function->GetParameter(4));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "([0]+[1]*(x-3.1)+[2]*(x-3.1)*(x-3.1))", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(5));
                function_BG->SetParameter(1, function->GetParameter(6));
                function_BG->SetParameter(2, function->GetParameter(7));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
        }

         void Crystall_Ball_fit_exp(TString options, TString name){

                int nb_param = 7;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                cout << "name BG func" << endl;
                cout << name_bg << endl;

                
                function = new TF1(name_function,"[0]*crystalball_function(x, [1], [2], [3], [4]) +  exp([5]+[6]*x)",min_fit, max_fit);
                //function = new TF1(name_function,"[0]*ROOT::Math::crystalball_function(x,[3], [4], [2], [1]) + ([5]+[6]*(x-3.1) + [7]*(x-3.1)*(x-3.1))",min_fit, max_fit);
                
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100.);
                function->SetParameter (1,3.096);
                function->SetParLimits (1,3.02, 3.2);
                function->SetParameter (2,0.04);
                function->SetParLimits (2,0.025, 0.15);
                function->SetParameter(3,0.75);
                function->SetParLimits(3,0.05,1.5);
                function->SetParameter(4,150);
                function->SetParLimits (4,100,500);
                //BG
                function->SetParameter(5, 7.);
                function->SetParLimits(5, 0.00, 100000.);
                function->SetParameter(6, -2.);
                function->SetParLimits(6, -100000., 0.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();

                function_Signal = new TF1(name_signal,"[0]*crystalball_function(x, [1], [2], [3], [4])", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetParameter(3, function->GetParameter(3));
                function_Signal->SetParameter(4, function->GetParameter(4));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "exp([0]+[1]*x)", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(5));
                function_BG->SetParameter(1, function->GetParameter(6));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
        }

        void Single_Gaussian_fit(TString options, TString name)
        {
                int nb_param = 5;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                cout << "name BG func" << endl;
                cout << name_bg << endl;

                function = new TF1(name_function, "([0])*exp(-0.5*pow((x-[1])/[2],2)) + exp([3]+[4]*x)", min_fit, max_fit);
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, init_amp_fit * 0.1, init_amp_fit * 100.);
                function->SetParameter(1, 3.096);
                function->SetParLimits(1, 3.02, 3.2);
                function->SetParameter(2, 0.04);
                function->SetParLimits(2, 0.025, 0.15);

                function->SetParameter(3, 7.);
                function->SetParLimits(3, 0.00, 100000.);
                function->SetParameter(4, -2.);
                function->SetParLimits(4, -100000., 0.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();

                function_Signal = new TF1(name_signal, "([0])*exp(-0.5*pow((x-[1])/[2],2)) ", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "exp([0]+[1]*x)", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(3));
                function_BG->SetParameter(1, function->GetParameter(4));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
                
        }

        void Single_Gaussian_Int_fit(TString options, TString name)
        {
                int nb_param = 5;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                cout << "name BG func" << endl;
                cout << name_bg << endl;

                function = new TF1(name_function, "([0]/(2.5066*[2]))*exp(-0.5*pow((x-[1])/[2],2)) + exp([3]+[4]*x)", min_fit, max_fit);
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, 0.000, 100000.);
                function->SetParameter(1, 3.096);
                function->SetParLimits(1, 3.02, 3.2);
                function->SetParameter(2, 0.04);
                function->SetParLimits(2, 0.025, 0.10);

                function->SetParameter(3, 7.);
                function->SetParLimits(3, 0.00, 100000.);
                function->SetParameter(4, -2.);
                function->SetParLimits(4, -100000., 0.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();

                function_Signal = new TF1(name_signal, "([0]/(2.5066*[2]))*exp(-0.5*pow((x-[1])/[2],2)) ", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "exp([0]+[1]*x)", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(3));
                function_BG->SetParameter(1, function->GetParameter(4));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
                
        }

        void Single_Gaussian_Int_fit_Pol_BG(TString options, TString name)
        {
                int nb_param = 6;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                cout << "name BG func" << endl;
                cout << name_bg << endl;

                function = new TF1(name_function, "([0]/(2.5066*[2]))*exp(-0.5*pow((x-[1])/[2],2)) + ([3]+[4]*(x-3.1)+[5]*(x-3.1)*(x-3.1))", min_fit, max_fit);
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, 0.000, 100000.);
                function->SetParameter(1, 3.096);
                function->SetParLimits(1, 3.02, 3.2);
                function->SetParameter(2, 0.04);
                function->SetParLimits(2, 0.025, 0.10);

                function->SetParameter(3, 7.);
                function->SetParLimits(3, 0.00, 100000.);
                function->SetParameter(4, -2.);
                //function->SetParLimits(4, -100., 100);
                function->SetParameter(5, 4.);
                function->SetParLimits(5, 0., 100000.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();

                function_Signal = new TF1(name_signal, "([0]/(2.5066*[2]))*exp(-0.5*pow((x-[1])/[2],2)) ", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "([0]+[1]*(x-3.1)+[2]*(x-3.1)*(x-3.1))", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(3));
                function_BG->SetParameter(1, function->GetParameter(4));
                function_BG->SetParameter(2, function->GetParameter(5));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
                
        }

        void Single_Gaussian_Int_fit_Pol_BG_V2(TString options, TString name)
        {
                int nb_param = 6;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                cout << "name BG func" << endl;
                cout << name_bg << endl;

                function = new TF1(name_function, "([0]/(2.5066*[2]))*exp(-0.5*pow((x-[1])/[2],2)) + ([3]+[5]*(x-[4])*(x-[4]))", min_fit, max_fit);
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, 0.000, 100000.);
                function->SetParameter(1, 3.096);
                function->SetParLimits(1, 3.02, 3.2);
                function->SetParameter(2, 0.04);
                function->SetParLimits(2, 0.025, 0.10);

                function->SetParameter(3, 7.);
                function->SetParLimits(3, 0.00, 100000.);
                function->SetParameter(4, 3.2);
                function->SetParLimits(4, 3.00, 3.5);
                function->SetParameter(5, 4.);
                function->SetParLimits(5, 0., 100000.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();

                function_Signal = new TF1(name_signal, "([0]/(2.5066*[2]))*exp(-0.5*pow((x-[1])/[2],2)) ", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "([0]+[2]*(x-[1])*(x-[1]))", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(3));
                function_BG->SetParameter(1, function->GetParameter(4));
                function_BG->SetParameter(2, function->GetParameter(5));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
                
        }

        void Double_Gaussian_Fit(TString options, TString name)
        {

                int nb_param = 7;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                function = new TF1(name_function, "([0])*exp(-0.5*pow((x-[1])/[2],2)) + ([3])*exp(-0.5*pow((x-[1])/[4],2)) + exp([5]+[6]*x)", min_fit, max_fit);
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, init_amp_fit * 0.2, init_amp_fit * 10.);
                function->SetParameter(1, 3.096);
                function->SetParLimits(1, 3.02, 3.2);
                function->SetParameter(2, 0.04);
                function->SetParLimits(2, 0.025, 0.15);

                function->SetParameter(3, init_amp_fit * 2);
                function->SetParLimits(3, 0.0, init_amp_fit * 5.);
                function->SetParameter(4, 0.03);
                function->SetParLimits(4, 0.020, 0.4);

                function->SetParameter(5, 7.);
                function->SetParLimits(5, 0.00, 100000.);
                function->SetParameter(6, -2.);
                function->SetParLimits(6, -100000., 0.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();

                function_Signal = new TF1(name_signal, "([0])*exp(-0.5*pow((x-[1])/[2],2)) + ([3])*exp(-0.5*pow((x-[1])/[4],2))", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetParameter(3, function->GetParameter(3));
                function_Signal->SetParameter(4, function->GetParameter(4));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "exp([0]+[1]*x)", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(3));
                function_BG->SetParameter(1, function->GetParameter(4));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
        }

        void Single_Gaussian_Fit_Flat_BG(TString options, TString name)
        {
                int nb_param = 5;
                cov_matrix.ResizeTo(nb_param, nb_param);

                TString name_function = name + "_fit_func";
                TString name_signal = name + "_sig_func";
                TString name_bg = name + "_bg_func";

                function = new TF1(name_function, "([0])*exp(-0.5*pow((x-[1])/[2],2)) + [3]+[4]*x", min_fit, max_fit);
                int init_amp_fit = (input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) > 0.0) ? input_Data_hist->GetBinContent(input_Data_hist->FindBin(3.096)) : 5;
                function->SetParameter(0, init_amp_fit);
                function->SetParLimits(0, init_amp_fit * 0.2, init_amp_fit * 10.);
                function->SetParameter(1, 3.096);
                function->SetParLimits(1, 3.02, 3.2);
                function->SetParameter(2, 0.04);
                function->SetParLimits(2, 0.025, 0.15);

                function->SetParameter(3, 7.);
                function->SetParLimits(3, 0.00, 100000.);
                function->SetParameter(4, -2.);
                function->SetParLimits(4, -100000., 0.0);

                input_Data_hist->Draw("e");
                TFitResultPtr fitResult = input_Data_hist->Fit(name_function, options);

                cov_matrix = fitResult->GetCovarianceMatrix();
                chi2 = fitResult->Chi2();
                NDF = fitResult->Ndf();
                // TMatrixDSym covarianceMatrix(function->GetNpar());
                // covarianceMatrix = fitResult->GetCovarianceMatrix();
                // cov_matrix = covarianceMatrix;

                function_Signal = new TF1(name_signal, "([0])*exp(-0.5*pow((x-[1])/[2],2)) ", min_fit, max_fit);
                function_Signal->SetParameter(0, function->GetParameter(0));
                function_Signal->SetParameter(1, function->GetParameter(1));
                function_Signal->SetParameter(2, function->GetParameter(2));
                function_Signal->SetLineColor(kGreen);
                function_Signal->SetFillColorAlpha(kGreen, 0.8);  // Set the transparency for the band
                function_Signal->SetLineWidth(3);

                function_BG = new TF1(name_bg, "exp([0]+[1]*x)", min_fit, max_fit);
                function_BG->SetParameter(0, function->GetParameter(3));
                function_BG->SetParameter(1, function->GetParameter(4));
                function_BG->SetLineColor(kBlue);
                function_BG->SetFillColorAlpha(kBlue, 0.8);
        }
};

#endif
