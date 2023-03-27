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

int ican2 = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);
}
void analysis() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   TFile * fo = new TFile( "dataPlots.root", "RECREATE" );

   auto * mVertexZ = new TH1F("mVertexZ", "VertexZ", 500, -100, 100);
   auto * mDeltaVertexZ = new TH1F("mDeltaVertexZ", "DeltaVertex", 500, 900, 1100);
   auto * mGRefMult = new TH1F("mGRefMult", "GRefMult", 500, -3, 3);
   auto * mZDCEast = new TH1F("mZDCEast", "ZDCEast", 500, 0, 500);
   auto * mZDCWest = new TH1F("mSDCWest", "ZDCWest", 500, 0, 500);

   auto * mChargeSum = new TH1F("ChargeSum", "Charge Sum", 500, -3, 3);
   auto * mPt = new TH1F("mPt", "Parent Transverse Momentum", 500, 0, 3);
   auto * mPhi = new TH1F("mPhi", "Parent Polar Angle (rad)", 500, -3.15, 3.15);
   auto * mEta = new TH1F("mEta", "Parent Pseudorapidity", 500, -6, 6);
   auto * mMass = new TH1F("mMass", "Parent Mass (GeV)", 500, 0.65, 0.9);
   auto * mRapidity = new TH1F("mRapidity", "Rapidity", 500, -2, 2);

    // Open the file containing the tree.
    TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2 ) + pow( pair->d2_mNSigmaPion, 2 );
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        Float_t mVertexZVal = pair->mVertexZ;
        Float_t mDeltaVertexZVal = pair->mDeltaVertexZ;        // delta vz
        UShort_t mGRefMultVal = pair->mGRefMult;              // global RefMult
        UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pair->mZDCWest;

        //Char_t    mChargeSumVal = pair->mChargeSum;
        Float_t   mRapidityVal = pair->mRapidity;
        

        if ( chipipi < 10 && dca1 < 1 && dca2 < 1){

            mVertexZ->Fill( mVertexZVal );
            mDeltaVertexZ->Fill( mDeltaVertexZVal );
            mGRefMult->Fill( mGRefMultVal );
            mZDCEast->Fill( mZDCEastVal );
            mZDCWest->Fill( mZDCWestVal );
            
            lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
            lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 ); 

            lv = lv1 + lv2;
            lvn = lv1 - lv2;
            int qs = pair->mChargeSum;


            if ( lv.M() > 0.65 && lv.M() < 0.9 ){
                mPt->Fill( lv.Pt() );
                mMass->Fill( lv.M() );
                mEta->Fill( lv.Eta() );
                mPhi->Fill( lv.Phi() );
                mChargeSum->Fill( qs );
                mRapidity->Fill( mRapidityVal );
                //ht->Fill( lv.Pt()*lv.Pt() );
            }

        }
    }

fo -> cd();

makeCanvas();
mVertexZ->SetLineColor(kBlack);
mVertexZ->Draw();
gPad->Print( "plot_mVertexZ.pdf" );

makeCanvas();
mDeltaVertexZ->SetLineColor(kBlack);
mDeltaVertexZ->Draw();
gPad->Print( "plot_mDeltaVertexZ.pdf" );

makeCanvas();
mGRefMult->SetLineColor(kBlack);
mGRefMult->Draw();
gPad->Print( "plot_mGRefMult.pdf" );

makeCanvas();
mZDCEast->SetLineColor(kBlack);
mZDCEast->Draw();
gPad->Print( "plot_mZDCEast.pdf" );

makeCanvas();
mZDCWest->SetLineColor(kBlack);
mZDCWest->Draw();
gPad->Print( "plot_mZDCWest.pdf" );

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();
gPad->Print( "plot_mMass.pdf" );

makeCanvas();
mPt->SetLineColor(kBlack);
mPt->Draw();
gPad->Print( "plot_mPt.pdf" );

makeCanvas();
mEta->SetLineColor(kBlack);
mEta->Draw();
gPad->Print( "plot_mEta.pdf" );

makeCanvas();
mPhi->SetLineColor(kBlack);
mPhi->Draw();
gPad->Print( "plot_mPhi.pdf" );

makeCanvas();
mRapidity->SetLineColor(kBlack);
mRapidity->Draw();
gPad->Print( "plot_mRapidity.pdf" );

makeCanvas();
mChargeSum->SetLineColor(kBlack);
mChargeSum->Draw();
gPad->Print( "plot_mChargeSum.pdf" );

fo->Write();
}



