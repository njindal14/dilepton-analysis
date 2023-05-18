//ttree reading and making some plots

#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

//#ifndef FEMTO_PAIR_H
//#define FEMTO_PAIR_H

int ican5 = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican5++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);
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



void eventvariables() {
    // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   TFile * fo = new TFile( "dataPlots.root", "RECREATE" );

   auto * mVertexZ = new TH1F("mVertexZ", "VertexZ", 500, -200, 200);
   auto * mDeltaVertexZ = new TH1F("mDeltaVertexZ", "#Delta VertexZ", 500, -1000, 1000);
   auto * mGRefMult = new TH1F("mGRefMult", "GRefMult", 700, 0, 700);
   auto * mZDCEast = new TH1F("mZDCEast", "ZDCEast", 1200, 0, 1200);
   auto * mZDCWest = new TH1F("mZDCWest", "ZDCWest", 1200, 0, 1200);
   //TF1 *peak1 = new TF1("peak1", "gaus", 25, 80);
   //TF1 *peak2 = new TF1("peak2", "gaus", 80, )

   //create TF1 to fit to the ZDC East
   auto *ZDCfunc = new TF1("fit", ZDCfitfunc, 30, 1200.0, 11);
   ZDCfunc->SetParameters(100000.0, 50.0, 15.0, 40000.0, 110.0, 30.0, 15000.0, 200.0, 30.0, 0.1, 5000.0);
   ZDCfunc->SetParNames("A1", "#lambda1", "#sigma1", "A2","#lambda2","#sigma2", "A3", "#lambda3", "#sigma3", "k", "p0");

   //set some parameter ranges
   ZDCfunc->SetParLimits(1,30,70);
   ZDCfunc->SetParLimits(2,0,30);
   ZDCfunc->SetParLimits(4,90,130);
   ZDCfunc->SetParLimits(5,25,35);
   ZDCfunc->SetParLimits(7,150,250);
   ZDCfunc->SetParLimits(8,20,50);
   ZDCfunc->SetParLimits(9,0,0.2);
   //smooth it out
   ZDCfunc->SetNpx(1000);
   ZDCfunc->SetLineWidth(4);



   

    TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    //TLorentzVector lv1, lv2, lv, lvn;

    while (myReader.Next()) {
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        Float_t mVertexZVal = pair->mVertexZ;
        Float_t mDeltaVertexZVal = pair->mDeltaVertexZ;        // delta vz
        UShort_t mGRefMultVal = pair->mGRefMult;              // global RefMult
        UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pair->mZDCWest;

        //Char_t    mChargeSumVal = pair->mChargeSum;
        //Float_t   mRapidityVal = pair->mRapidity;


        //if ( chiee < 10 && dca1 < 1 && dca2 < 1){

            mVertexZ->Fill( mVertexZVal );
            mDeltaVertexZ->Fill( mDeltaVertexZVal );
            mGRefMult->Fill( mGRefMultVal );
            if(mZDCEastVal > 30){
                mZDCEast->Fill( mZDCEastVal );
            }
            if(mZDCWestVal > 30){
                mZDCWest->Fill( mZDCWestVal );
            }

    //}
}
    
fo -> cd();

/*
makeCanvas();
mVertexZ->SetLineColor(kBlack);
gPad->SetLogy();
mVertexZ->GetXaxis()->SetTitle("Z Vertex Location (cm)");
mVertexZ->GetYaxis()->SetTitle("Counts");
mVertexZ->Draw();
gPad->Print( "plot_mVertexZ.png" );

makeCanvas();
mDeltaVertexZ->SetLineColor(kRed);
gPad->SetLogy();
mDeltaVertexZ->GetXaxis()->SetTitle("#Delta Z Vertex Location (cm)");
mDeltaVertexZ->GetYaxis()->SetTitle("Counts");
mDeltaVertexZ->Draw();
gPad->Print( "plot_mDeltaVertexZ.png" );

makeCanvas();
mGRefMult->SetLineColor(kBlack);
gPad->SetLogy();
mGRefMult->GetXaxis()->SetTitle("GRefMult");
mGRefMult->GetYaxis()->SetTitle("Counts");
mGRefMult->Draw();
gPad->Print( "plot_mGRefMult.png" );
*/

makeCanvas();
mZDCEast->SetLineColor(kBlack);
mZDCEast->GetXaxis()->SetTitle("ZDC East");
mZDCEast->GetYaxis()->SetTitle("Counts");
mZDCEast->Draw();
mZDCEast->Fit("fit", "", "", 30,1200);
//gStyle->SetLineWidth(4);
gStyle->SetOptFit(1111);

//get parameters to later plot each one individually
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
cout << A1;
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
gPad->Print( "plot_mZDCEast1Fit.png" );

makeCanvas();
mZDCWest->SetLineColor(kRed);
mZDCWest->SetLineWidth(5);
mZDCWest->GetXaxis()->SetTitle("ZDC West");
mZDCWest->GetYaxis()->SetTitle("Counts");
mZDCWest->Draw();
gPad->Print( "plot_mZDCWest.png" );
}

