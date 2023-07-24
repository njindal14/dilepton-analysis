#ifndef FEMTO_PAIR_H
#define FEMTO_PAIR_H

#include "TObject.h"

int ican = 0;
void makeCan() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);
}
void rhopAu() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   //TFile * fo = new TFile( "data.root", "RECREATE" );

   auto * mVertexZ = new TH1F("mVertexZ", "", 500, 0, 3);
   auto * mDeltaVertexZ = new TH1F("mDeltaVertexZ", "", 500, 0, 5);
   auto * mGRefMult = new TH1F("mGRefMult", "", 500, 0, 5);
   auto * mZDCEast = new TH1F("mZDCEast", "", 500, 0, 1);
   auto * mZDCWest = new TH1F("mSDCWest", "", 500, 0, 1);

   auto * mChargeSum = new TH1F("ChargeSum", "", 500, 0, 1);
   auto * mPt = new TH1F("mPt", "", 500, 0, 3);
   auto * mEta = new TH1F("mEta", "", 500, 0, 6);
   auto * mMass = new TH1F("mMass", "", 500, 0, 5);
   auto * mRapidity = new TH1F("mRapidity", "", 500, 0, 10);

    // Open the file containing the tree.
    TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2 ) + pow( pair->d1_mNSigmaPion, 2 );
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        if ( chipipi < 10 && dca1 < 1 && dca2 < 1){
            
            lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
            lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 ); 

            lv = lv1 + lv2;
            lvn = lv1 - lv2;
            int qs = pair->mChargeSum;

            hmass->Fill( lv.M() );

            if ( lv.M() > 0.65 && lv.M() < 0.9 ){
                hpt->Fill( lv.Pt() );
                ht->Fill( lv.Pt()*lv.Pt() );
            }

        }
    }

//fo -> cd();
makeCan();
hulsmass->SetLineColor(kBlack);
hulsmass->Draw();
}



