#include "TLegend.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"

//for plotting
int ican10 = 0;
void makeCanvas5() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican10++ ), "", 900, 600 );
    can->SetTopMargin(0.1);
    can->SetBottomMargin(0.1);
    can->SetRightMargin(0.05);

}

static const double hbarc    = 0.197327053;
static const double Znu       = 79;
static const double RNuc       = 6.38;
static const double Anu       = 197;
double A=Anu;
double c = 3;


double formFactor(double * x, double *par)
{
	//par[0] is R, x is t (pt^2)
    //par[0] = pow(A, 1./3.) * 1.16 * (1. - 1.16 * pow(A, -2. / 3.));
	double q    = sqrt(x[0]/2)*c;
	double arg1 = q * par[0] / hbarc;
	double arg2 = hbarc / (q * 1.16 * (1. - 1.16 * pow(A, -2. / 3.)));
	double sph  = (sin(arg1) - arg1 * cos(arg1)) * 3. * arg2 * arg2 * arg2 / double(A);
	double a0   = 0.70;  // [fm]
	return sph / (1. + (a0 * a0 * x[0]/sqrt(2)*c) / (hbarc * hbarc));
}

void formfactor() {
    auto * formfactorLowR = new TF1("ffLowR", formFactor, 0, .02, 1);
    auto * formfactorMidR = new TF1("ffMidR", formFactor, 0, .02, 1);
    auto * formfactorHighR = new TF1("ffHighR", formFactor, 0, .02, 1);

    formfactorLowR->SetParameters(3);
    formfactorMidR->SetParameters(5);
    formfactorHighR->SetParameters(7);

    makeCanvas5();
    formfactorHighR->SetTitle("Form Factors for Various Radii (fm)");
    formfactorHighR->GetXaxis()->SetTitle("b [fm]");
    formfactorHighR->GetYaxis()->SetTitle("#rho(b) [fm^{-2}]");
    formfactorHighR->SetLineColor(kBlack);
    formfactorHighR->Draw();
    formfactorMidR->SetLineColor(kRed);
    formfactorMidR->Draw("same");
    formfactorLowR->SetLineColor(kGreen);
    formfactorLowR->Draw("same");
    auto * ffLegend = new TLegend(0.79,0.3,.99,0.55);
    ffLegend->AddEntry(formfactorHighR, "R = 7 fm");
    ffLegend->AddEntry(formfactorMidR, "R = 5 fm");
    ffLegend->AddEntry(formfactorLowR, "R = 3 fm");
    ffLegend->SetTextSize(0.03);
    ffLegend->Draw("same");
    gPad->Print("plot_formfactors.png");

}


