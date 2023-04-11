//ttree reading and making some plots

#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"

//#ifndef FEMTO_PAIR_H
//#define FEMTO_PAIR_H

int ican5 = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican5++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);
}
void eventvariables() {
    // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   TFile * fo = new TFile( "dataPlots.root", "RECREATE" );

   auto * mVertexZ = new TH1F("mVertexZ", "VertexZ", 500, -100, 100);
   auto * mDeltaVertexZ = new TH1F("mDeltaVertexZ", "#Delta VertexZ", 500, 900, 1100);
   auto * mGRefMult = new TH1F("mGRefMult", "GRefMult", 500, -3, 3);
   auto * mZDCEast = new TH1F("mZDCEast", "ZDCEast", 750, 0, 750);
   auto * mZDCWest = new TH1F("mZDCWest", "ZDCWest", 750, 0, 750);


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


        if ( chiee < 10 && dca1 < 1 && dca2 < 1){

            mVertexZ->Fill( mVertexZVal );
            mDeltaVertexZ->Fill( mDeltaVertexZVal );
            mGRefMult->Fill( mGRefMultVal );
            mZDCEast->Fill( mZDCEastVal );
            mZDCWest->Fill( mZDCWestVal );

    }
}
    
fo -> cd();

makeCanvas();
mVertexZ->SetLineColor(kBlack);
mVertexZ->GetXaxis()->SetTitle("Z Vertex Location (cm)");
mVertexZ->GetYaxis()->SetTitle("Counts");
mVertexZ->Draw();
gPad->Print( "plot_mVertexZ.png" );

makeCanvas();
mDeltaVertexZ->SetLineColor(kRed);
mDeltaVertexZ->GetXaxis()->SetTitle("#Delta Z Vertex Location (cm)");
mDeltaVertexZ->GetYaxis()->SetTitle("Counts");
mDeltaVertexZ->Draw();
gPad->Print( "plot_mDeltaVertexZ.png" );

makeCanvas();
mGRefMult->SetLineColor(kBlack);
mGRefMult->GetXaxis()->SetTitle("GRefMult");
mGRefMult->GetYaxis()->SetTitle("Counts");
mGRefMult->Draw();
gPad->Print( "plot_mGRefMult.png" );

makeCanvas();
mZDCEast->SetLineColor(kBlack);
mZDCEast->GetXaxis()->SetTitle("ZDC East");
mZDCEast->GetYaxis()->SetTitle("Counts");
mZDCEast->Draw();
gPad->Print( "plot_mZDCEast.png" );

makeCanvas();
mZDCWest->SetLineColor(kRed);
mZDCWest->GetXaxis()->SetTitle("ZDC West");
mZDCWest->GetYaxis()->SetTitle("Counts");
mZDCWest->Draw();
gPad->Print( "plot_mZDCWest.png" );
}

