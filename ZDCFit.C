
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"
#include <cmath>

//#ifndef FEMTO_PAIR_H
//#define FEMTO_PAIR_H

int ican6 = 0;
void makeCanvas3() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican6++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.3);
}

double ZDCfitfunc(double *x, double *par){
    double fitval;
    if((x[0]-par[7])/par[8] > par[9]){
        fitval = par[0]/(par[2]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + 
            par[3]/(par[5]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[4])/par[5], 2)) +
                par[6]/(par[8]*sqrt(2*M_PI)) *exp(pow(par[9], 2)/2 - par[9]*(x[0]-par[7])/par[8]) + par[10];                   
    }
    else {
        fitval = par[0]/(par[2]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + 
            par[3]/(par[5]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[4])/par[5], 2)) +
                par[6]/(par[8]*sqrt(2*M_PI)) *exp(-0.5*pow((x[0]-par[7])/par[8], 2)) +par[10];
    }
    return fitval;
}

double ptfitfunc(double *x, double *par){
    double fitval;
    if((x[0] - par[1])/par[2] > par[3]){
        fitval = par[0]/(par[2]*sqrt(2*M_PI)) * exp(pow(par[3], 2)/2 - par[3]*(x[0]-par[1])/par[2]) + par[4];
    }
    else{
        fitval = par[0]/(par[2]*sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + par[4];
    }
    return fitval;
}

double calcNeutrons(double x, double mean1, double sigma1, double mean2, double sigma2, double mean3, double sigma3){
    
    double nNeutrons;
    double z1 = (x-mean1)/sigma1;
    double z2 = (x-mean2)/sigma2;
    double z3 = (x-mean3)/sigma3;
    double probability1 = 0.5 * (1 + std::erf(z1 / std::sqrt(2.0)));
    double probability2 = 0.5 * (1 + std::erf(z2 / std::sqrt(2.0)));
    double probability3 = 0.5 * (1 + std::erf(z3 / std::sqrt(2.0)));
    nNeutrons = probability1 + 2*probability2 + 3*probability3;
    return nNeutrons;
    
}

void ZDCFit() {
    // Create a histogram for the values we read.
    TH1F("h1", "ntuple", 100, -4, 4);

    auto * mZDCEast = new TH1F("mZDCEast, ee Events (PID + TOF + Mass Cut)", "ZDCEast", 400, 0, 1200);
    auto * mZDCWest = new TH1F("mZDCWest, ee Events (PID + TOF + Mass Cut)", "ZDCWest", 400, 0, 1200);

    auto * ZDCsSquared = new TH1F("ZDC's Squared", "ZDC's Squared", 500, 0, 2000);

    auto * mZDCEast1n1n = new TH1F("mZDCEast Separated, ee Events", "ZDCEast Separated", 1200, 0, 1200);
    auto * mZDCWest1n1n = new TH1F("mZDCWest Separated, ee Events", "ZDCWest Separated", 1200, 0, 1200);

    auto * mZDCEast2n2n = new TH1F("mZDCEast Separated, ee Events", "ZDCEast Separated", 1200, 0, 1200);
    auto * mZDCWest2n2n = new TH1F("mZDCWest Separated, ee Events", "ZDCWest Separated", 1200, 0, 1200);

    auto * mZDCEast3n3nPlus = new TH1F("mZDCEast3n3nPlus", "ZDCEast3n3nPlus", 1200, 0, 1200);
    auto * mZDCWest3n3nPlus = new TH1F("mZDCWest3n3nPlus", "ZDCWest3n3nPlus", 1200, 0, 1200);

    auto * mPt1n1nUnlikeSign = new TH1F("Pair Transverse Momentum, ee Events", "Pair Transverse Momentum", 100, 0, .15);
    auto * mPt2n2nUnlikeSign = new TH1F("Pair Transverse Momentum, ee Events", "Pair Transverse Momentum", 100, 0, .15);
    auto * mPt3n3nPlusLikeSign = new TH1F("Pair Transverse Momentum, ee Events", "Pair Transverse Momentum", 100, 0, .15);
    auto * mPt3n3nPlusUnlikeSign = new TH1F("Pair Transverse Momentum, ee Events", "Pair Transverse Momentum", 100, 0, .15);


    auto * numNeutronsEast = new TH1F("Continuous Number of Neutrons East", "Neutrons East", 200, 0, 5);
    auto * mPt = new TH1F("mPt, PID, M_{ee} + TOF Cuts", "Parent Transverse Momentum (Gev/c)", 500, 0, 1);

    auto * numNeutronsvPtEast = new TH2F("Neutrons East Vs. P_{T}", "Neutrons East Vs P_{T}", 500, 0, .15, 20, 0,5);

    auto * numNeutronsWest = new TH1F("Continuous Number of Neutrons West", "Neutrons West", 200, 0, 5);
    auto * numNeutronsvPtWest = new TH2F("Neutrons West Vs. P_{T}", "Neutrons West Vs P_{T}", 500, 0, .15, 20, 0,5);








    auto * ZDC2D = new TH2F("ZDC Heat Map", "ZDC Heat Map", 1200, 0, 1200, 1200, 0, 1200);


    auto *ZDCfunc = new TF1("fit", ZDCfitfunc, 30, 1200.0, 11);
    ZDCfunc->SetParameters(100000.0, 50.0, 15.0, 40000.0, 110.0, 30.0, 15000.0, 200.0, 30.0, 0.1, 5000.0);
    ZDCfunc->SetParNames("A1", "#lambda1", "#sigma1", "A2","#lambda2","#sigma2", "A3", "#lambda3", "#sigma3", "k", "p0");

    auto *ptfunc = new TF1("ptfit", ptfitfunc, 0,0.1, 5);
    auto *ptfunc2 = new TF1("ptfit2", ptfitfunc, 0,0.1,5);
    ptfunc->SetParameters(2, .025, .01, .01, .01);
    ptfunc->SetParNames("a","mean", "#sigma","k", "p0");
    //ptfunc->SetParLimits(1,.02,.03);
    
    //ptfunc->SetParLimits(3, -1, 1);

    ptfunc2->SetParameters(1, .025, .01, .01, .01);
    ptfunc2->SetParNames("a2","mean2", "#sigma2","k2", "p02");
    //ptfunc2->SetParLimits(1,.02,.03);

    //set some parameter ranges
   ZDCfunc->SetParLimits(1,30,70);
   ZDCfunc->SetParLimits(2,10,30);
   ZDCfunc->SetParLimits(4,90,130);
   ZDCfunc->SetParLimits(5,25,35);
   ZDCfunc->SetParLimits(7,150,250);
   ZDCfunc->SetParLimits(8,20,50);
   ZDCfunc->SetParLimits(9,0,0.2);
   //smooth it out
   ZDCfunc->SetNpx(1000);
   ZDCfunc->SetLineWidth(4);

   ptfunc->SetNpx(1000);
   ptfunc->SetLineWidth(4);

   ptfunc2->SetNpx(1000);
   ptfunc2->SetLineWidth(4);


    TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;


    //TLorentzVector lv1, lv2, lv, lvn;
    int counter1 = 0;
    int counter2 = 0;
    int counter3 = 0;

    while (myReader.Next()) {
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pair->mZDCWest;
        UShort_t mZDCEastVal2 = pair->mZDCEast;
        UShort_t mZDCWestVal2 = pair->mZDCWest;
        Float_t mPtVal = pair->mPt;
        int chargesumval = pair->mChargeSum;

        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        //get parent total momentum for each track
        Float_t p1 = lv1.P();
        Float_t p2 = lv2.P();
        //square it
        Float_t p1_2 = pow(p1,2);
        Float_t p2_2 = pow(p2,2);
        Float_t mVertexZVal = pair->mVertexZ;
        UShort_t mGRefMultVal = pair->mGRefMult;  
        Float_t Tof1 = pair->d1_mTof;     
        Float_t Tof2 = pair->d2_mTof;  
        Float_t len1 = pair->d1_mLength;
        Float_t len2 = pair->d2_mLength; 
        //Float_t mPtVal = pair->mPt;
        Float_t dTofVal = Tof1 - Tof2;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp1 - texp2;
        Float_t ddTofVal = dTofVal - dTofexpVal;
        


        //PID and TOF cuts
        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
       pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {

            if(chiee < 10 && 3*chiee < chipipi && ddTofVal < 0.5 && ddTofVal > -0.5){


                //get the right phase space for BW process, cut on mass of parent and calculate
                //neutrons as a continous variable from zdc fit numbers
                if(mZDCWestVal > 30 && mZDCEastVal > 30 ){
                    if(lv.M() < 0.8 && lv.M() > 0.5){
                        mPt->Fill( mPtVal ); 
                        mZDCWest->Fill( mZDCWestVal );
                        double numNeutronswest = calcNeutrons(mZDCWestVal, 54.12, 13, 105.1, 25.3, 150.8, 33);
                        numNeutronsWest->Fill(numNeutronswest);
                        numNeutronsvPtWest->Fill(mPtVal, numNeutronswest);

                        mZDCEast->Fill( mZDCEastVal );
                        double numNeutronseast = calcNeutrons(mZDCEastVal, 52.12, 12.75, 109.4, 25.3, 172.8, 20);
                        numNeutronsEast->Fill(numNeutronseast);
                        numNeutronsvPtEast->Fill(mPtVal, numNeutronseast);
                    }
                }

                ZDCsSquared->Fill(sqrt(pow(mZDCEastVal, 2) + pow(mZDCWestVal, 2)));
                //1n1n cuts (hard numbers, just to check what it looks like)
                if(mZDCEastVal > 30 && mZDCWestVal > 30 && mZDCEastVal < 80
                    && mZDCWestVal < 80){
                        mZDCEast1n1n->Fill( mZDCEastVal);          
                        mZDCWest1n1n->Fill( mZDCWestVal);
                        counter1 +=1;
                        mPt1n1nUnlikeSign->Fill(mPtVal);
                    }

                //2n2n cuts
                if(mZDCEastVal > 30 && mZDCWestVal > 30 && mZDCEastVal < 155 && mZDCEastVal > 80
                    && mZDCWestVal < 155 && mZDCWestVal > 80){
                        mZDCWest2n2n->Fill(mZDCWestVal);
                        mZDCEast2n2n->Fill(mZDCEastVal);
                        counter2 +=1;
                        mPt2n2nUnlikeSign->Fill(mPtVal);
                    }

                //3n3nPlus cuts
                if(mZDCEastVal > 155 && mZDCWestVal > 155){
                        mZDCWest3n3nPlus->Fill(mZDCWestVal);
                        mZDCEast3n3nPlus->Fill(mZDCEastVal);
                        counter3 +=1;
                        mPt3n3nPlusUnlikeSign->Fill(mPtVal);
                    }
                    ZDC2D->Fill(mZDCEastVal2, mZDCWestVal2);
            
                
            }
        }


    }
    cout << counter1;
    cout << "first";
    cout << counter2;
    cout << "second";
    cout << counter3;

    //projections to check my cuts
    auto * projEast = ZDC2D->ProjectionX("ZDCEastcut1nw", 30, 90);
    auto * projWest = ZDC2D->ProjectionY("ZDCWestcut1ne", 30, 90);


makeCanvas3();
projEast->SetLineColor(kBlack);
projEast->Fit("fit","","",30,1200);
gStyle->SetOptFit(1111);
projEast->Draw();
gPad->Print( "plots/plot_mZDC2DEastProjection1n.png" );


makeCanvas3();
projWest->SetLineColor(kRed);
projWest->Fit("fit","","",30,1200);
gStyle->SetOptFit(1111);
projWest->Draw();
gPad->Print( "plots/plot_mZDC2DWestProjection1n.png" );


//plot and fit ZDC's
makeCanvas3();
mZDCEast->SetLineColor(kBlack);
mZDCEast->GetXaxis()->SetTitle("ZDC East");
mZDCEast->GetYaxis()->SetTitle("Counts");
mZDCEast->Draw();
mZDCEast->Fit("fit", "", "", 30,1200);
gStyle->SetOptFit(1111);
//get parameters to later plot each one individually for ZDCEast
double A1 = ZDCfunc->GetParameter(0);
double A1e = ZDCfunc->GetParError(0);
double l1 = ZDCfunc->GetParameter(1);
double l1e = ZDCfunc->GetParError(1);
double s1 = ZDCfunc->GetParameter(2);
double s1e = ZDCfunc->GetParError(2);
double A2 = ZDCfunc->GetParameter(3);
double A2e = ZDCfunc->GetParError(3);
double l2 = ZDCfunc->GetParameter(4);
double l2e = ZDCfunc->GetParError(4);
double s2 = ZDCfunc->GetParameter(5);
double s2e = ZDCfunc->GetParError(5);
double A3 = ZDCfunc->GetParameter(6);
double A3e = ZDCfunc->GetParError(6);
double l3 = ZDCfunc->GetParameter(7);
double l3e = ZDCfunc->GetParError(7);
double s3 = ZDCfunc->GetParameter(8);
double s3e = ZDCfunc->GetParError(8);
double k3 = ZDCfunc->GetParameter(9);
double k3e = ZDCfunc->GetParError(9);
double p0 = ZDCfunc->GetParameter(10);
double p0e = ZDCfunc->GetParError(10);
auto *g1 = new TF1("g1", "[0]*exp(-0.5 * ((x-[1])/[2]) * ((x-[1])/[2]))", 0, 100);
g1->SetParameter(0,A1/(s1*sqrt(2*M_PI)));
//cout << A1;
g1->SetParameter(1,l1);
g1->SetParameter(2,s1);
auto *g2 = new TF1("g2", "[3]*exp(-0.5 * ((x-[4])/[5]) * ((x-[4])/[5]))", 0, 200);
g2->SetParameter(3,A2/(s2*sqrt(2*M_PI)));
g2->SetParameter(4,l2);
g2->SetParameter(5,s2);
double cutoff = s3*k3 + l3;
auto *g3 = new TF1("g3", "[6]* exp(-0.5 * ((x-[7])/[8]) * ((x-[7])/[8]) )", 0, cutoff+3);
g3->SetParameter(6, A3/(s3*sqrt(2*M_PI)));
g3->SetParameter(7, l3);
g3->SetParameter(8, s3);
auto *g4 = new TF1("g4", "[9]* exp([11]*[11]/2 - [11]*((x-[10])/[12]))", cutoff-3, 1200);
g4->SetParameter(9, A3/(s3*sqrt(2*M_PI)));
g4->SetParameter(10, l3);
g4->SetParameter(11, k3);
g4->SetParameter(12, s3);
auto *g5 = new TF1("g5", "[13]", 0, 1200);
g5->SetParameter(13, p0);
g1->SetLineColor(kBlue);
g1->SetLineWidth(4);
g1->Draw("same");
g2->SetLineColor(kBlack);
g2->SetLineWidth(4);
g2->Draw("same");
g3->SetLineColor(kPink);
g3->SetLineWidth(4);
g3->Draw("same");
g4->SetLineColor(kPink);
g4->SetLineWidth(4);
g4->Draw("same");
g5->SetLineColor(kGray);
g5->SetLineWidth(4);
g5->Draw("same");
gPad->Print( "plots/plot_mZDCEast2FitBW.png" );






makeCanvas3();
mZDCWest->SetLineColor(kBlack);
mZDCWest->GetXaxis()->SetTitle("ZDC West");
mZDCWest->GetYaxis()->SetTitle("Counts");
mZDCWest->Draw();
mZDCWest->Fit("fit", "", "", 30,1200);
gStyle->SetOptFit(1111);
//get parameters to later plot each one individually for ZDCWest
double A1w = ZDCfunc->GetParameter(0);
double A1we = ZDCfunc->GetParError(0);
double l1w = ZDCfunc->GetParameter(1);
double l1we = ZDCfunc->GetParError(1);
double s1w = ZDCfunc->GetParameter(2);
double s1we = ZDCfunc->GetParError(2);
double A2w = ZDCfunc->GetParameter(3);
double A2we = ZDCfunc->GetParError(3);
double l2w = ZDCfunc->GetParameter(4);
double l2we = ZDCfunc->GetParError(4);
double s2w = ZDCfunc->GetParameter(5);
double s2we = ZDCfunc->GetParError(5);
double A3w = ZDCfunc->GetParameter(6);
double A3we = ZDCfunc->GetParError(6);
double l3w = ZDCfunc->GetParameter(7);
double l3we = ZDCfunc->GetParError(7);
double s3w = ZDCfunc->GetParameter(8);
double s3we = ZDCfunc->GetParError(8);
double k3w = ZDCfunc->GetParameter(9);
double k3we = ZDCfunc->GetParError(9);
double p0w = ZDCfunc->GetParameter(10);
double p0we = ZDCfunc->GetParError(10);
auto *g1w = new TF1("g1w", "[0]*exp(-0.5 * ((x-[1])/[2]) * ((x-[1])/[2]))", 0, 100);
g1w->SetParameter(0,A1w/(s1w*sqrt(2*M_PI)));
//cout << A1;
g1w->SetParameter(1,l1w);
g1w->SetParameter(2,s1w);
auto *g2w = new TF1("g2w", "[3]*exp(-0.5 * ((x-[4])/[5]) * ((x-[4])/[5]))", 0, 200);
g2w->SetParameter(3,A2w/(s2w*sqrt(2*M_PI)));
g2w->SetParameter(4,l2w);
g2w->SetParameter(5,s2w);
double cutoffw = s3w*k3w + l3w;
auto *g3w = new TF1("g3w", "[6]* exp(-0.5 * ((x-[7])/[8]) * ((x-[7])/[8]) )", 0, cutoffw+3);
g3w->SetParameter(6, A3w/(s3w*sqrt(2*M_PI)));
g3w->SetParameter(7, l3w);
g3w->SetParameter(8, s3w);
auto *g4w = new TF1("g4w", "[9]* exp([11]*[11]/2 - [11]*((x-[10])/[12]))", cutoffw-3, 1200);
g4w->SetParameter(9, A3w/(s3w*sqrt(2*M_PI)));
g4w->SetParameter(10, l3w);
g4w->SetParameter(11, k3w);
g4w->SetParameter(12, s3w);
auto *g5w = new TF1("g5w", "[13]", 0, 1200);
g5w->SetParameter(13, p0w);
g1w->SetLineColor(kBlue);
g1w->SetLineWidth(4);
g1w->Draw("same");
g2w->SetLineColor(kBlack);
g2w->SetLineWidth(4);
g2w->Draw("same");
g3w->SetLineColor(kPink);
g3w->SetLineWidth(4);
g3w->Draw("same");
g4w->SetLineColor(kPink);
g4w->SetLineWidth(4);
g4w->Draw("same");
g5w->SetLineColor(kGray);
g5w->SetLineWidth(4);
g5w->Draw("same");
gPad->Print( "plots/plot_mZDCWest2FitBW.png" );



//color map of ZDC's
makeCanvas3();
//gPad->SetLogz();
ZDC2D->SetContour(100);
ZDC2D->GetXaxis()->SetTitle("ADC ZDC East");
ZDC2D->GetYaxis()->SetTitle("ADC ZDC West");
gStyle->SetPalette(1);
ZDC2D->Draw("colz");
gStyle->SetPalette(1);
gPad->Print( "plots/ZDC2DWide.png" );

makeCanvas3();
mZDCEast1n1n->SetLineColor(kBlack);
mZDCEast1n1n->GetXaxis()->SetTitle("ADC ZDC East");
mZDCEast1n1n->GetYaxis()->SetTitle("Counts");
mZDCEast1n1n->Draw();

mZDCEast2n2n->SetLineColor(kRed);
mZDCEast2n2n->Draw("same");

mZDCEast3n3nPlus->SetLineColor(kGreen);
mZDCEast3n3nPlus->Draw("same");
auto * legend = new TLegend(0.75,0.5,.95,0.6);
legend->SetHeader("Legend");
legend->AddEntry(mZDCEast1n1n,"1n1n Peak","l");
legend->AddEntry(mZDCEast2n2n,"2n2n Peak","l");
legend->AddEntry(mZDCEast3n3nPlus, "3n3nPlus Peak", "l");
legend->Draw("same");
gPad->Print( "plots/ZDCEastSeparated2VisualCutsBW.png" );



makeCanvas3();
mZDCWest1n1n->SetLineColor(kBlack);
mZDCWest1n1n->GetXaxis()->SetTitle("ADC ZDC West");
mZDCWest1n1n->GetYaxis()->SetTitle("Counts");
mZDCWest1n1n->Draw();

mZDCWest2n2n->SetLineColor(kRed);
mZDCWest2n2n->Draw("same");

mZDCWest3n3nPlus->SetLineColor(kGreen);
mZDCWest3n3nPlus->Draw("same");
auto * legend2 = new TLegend(0.75,0.5,.95,0.6);
legend2->SetHeader("Legend");
legend2->AddEntry(mZDCWest1n1n,"1n1n Peak","l");
legend2->AddEntry(mZDCWest2n2n,"2n2n Peak","l");
legend2->AddEntry(mZDCWest3n3nPlus, "3n3nPlus Peak", "l");
legend2->Draw("same");
gPad->Print( "plots/ZDCWestSeparatedVisualCutsBW.png" );


makeCanvas3();
mPt2n2nUnlikeSign->SetLineColor(kRed);
mPt2n2nUnlikeSign->GetXaxis()->SetTitle("Pt (GeV/c)");
mPt2n2nUnlikeSign->GetYaxis()->SetTitle("Counts (Normalized)");
mPt2n2nUnlikeSign->Scale(1/mPt2n2nUnlikeSign->GetEntries());
mPt2n2nUnlikeSign->Draw();

mPt3n3nPlusUnlikeSign->SetLineColor(kGreen);
mPt3n3nPlusUnlikeSign->GetXaxis()->SetTitle("Pt (GeV/c)");
mPt3n3nPlusUnlikeSign->GetYaxis()->SetTitle("Counts");
mPt3n3nPlusUnlikeSign->Scale(1/mPt3n3nPlusUnlikeSign->GetEntries());
mPt3n3nPlusUnlikeSign->Draw("same");
mPt3n3nPlusLikeSign->SetLineColor(kBlue);

TH1F * diff = (TH1F*)mPt3n3nPlusUnlikeSign->Clone();
diff->Add(mPt3n3nPlusLikeSign, -1);
diff->SetLineColor(kPink);

mPt1n1nUnlikeSign->SetLineColor(kBlack);
mPt1n1nUnlikeSign->Scale(1/mPt1n1nUnlikeSign->GetEntries());
mPt1n1nUnlikeSign->Draw("same");

auto * legend3 = new TLegend(0.75,0.3,.95,0.5);
legend3->SetHeader("Legend");
legend3->AddEntry(mPt1n1nUnlikeSign,"Pt 1n1n Peak Unlike Sign","l");
legend3->AddEntry(mPt2n2nUnlikeSign,"Pt 2n2n Peak Unlike Sign","l");
//legend3->AddEntry(mPt3n3nPlusLikeSign, "Pt 3n3nPlus Peak Like Sign", "l");
legend3->AddEntry(mPt3n3nPlusUnlikeSign, "Pt 3n3nPlus Peak Unlike Sign", "l");
//legend3->AddEntry(diff, "Pt 3n3nPlus Peak unlike-like", "l");

legend3->Draw("same");
//gPad->SetLogy();
gPad->Print( "plots/PtXnXnSeparatedVisualCutsbySignBWNormalized.png" );

makeCanvas3();
ZDCsSquared->SetLineColor(kBlack);
ZDCsSquared->GetXaxis()->SetTitle("sqrt(ZDCEast^2+ZDCWest^2)");
ZDCsSquared->GetYaxis()->SetTitle("Counts");
ZDCsSquared->Draw();
gPad->Print("plots/plot_ZDCsSquared.png");


makeCanvas3();
numNeutronsEast->SetLineColor(kBlack);
numNeutronsEast->GetXaxis()->SetTitle("Number of Neutrons ZDC East");
numNeutronsEast->GetYaxis()->SetTitle("Counts");
numNeutronsEast->Draw();
gPad->Print("plots/neutronsEastWeighted.png");

makeCanvas3();
//gPad->SetLogz();
numNeutronsvPtEast->SetContour(100);
numNeutronsvPtEast->GetXaxis()->SetTitle("P_{T}");
numNeutronsvPtEast->GetYaxis()->SetTitle("Number of Neutrons (Weighted)");
numNeutronsvPtEast->Draw("colz");
gStyle->SetPalette(1);
gPad->Print("plots/neutronseastvPt.png");


auto *numNeutronsEastvPtProjection = numNeutronsvPtEast->ProfileY("numNeutronsEastvPtProjection", 1, -1);
auto *numNeutronsWestvPtProjection = numNeutronsvPtWest->ProfileY("numNeutronsWestvPtProjection", 1, -1);


makeCanvas3();
numNeutronsEastvPtProjection->SetTitle("P_{T} Vs. Weighted East Neutron Mult");
numNeutronsEastvPtProjection->GetXaxis()->SetTitle("Neutron Multiplicity East");
numNeutronsEastvPtProjection->GetYaxis()->SetTitle("P_{T}");
numNeutronsEastvPtProjection->SetLineColor(kBlack);
numNeutronsEastvPtProjection->Draw();
gPad->Print( "plots/neutronmulteastvPt.png");

makeCanvas3();
numNeutronsWestvPtProjection->SetTitle("P_{T} Vs. Weighted West Neutron Mult");
numNeutronsWestvPtProjection->GetXaxis()->SetTitle("Neutron Multiplicity West");
numNeutronsWestvPtProjection->GetYaxis()->SetTitle("P_{T}");
numNeutronsWestvPtProjection->SetLineColor(kRed);
numNeutronsWestvPtProjection->Draw();
gPad->Print( "plots/neutronmultwestvPt.png");





}
