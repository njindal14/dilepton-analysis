
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

int ican3 = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican3++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);

}
void eeAnalysis() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);


   //TFile * fo = new TFile( "dataPlots.root", "RECREATE" );

   auto * mVertexZ = new TH1F("hist_vz", "VertexZ", 500, -100, 100);
   auto * mDeltaVertexZ = new TH1F("mDeltaVertexZ", "DeltaVertex", 500, 900, 1100);
   auto * mGRefMult = new TH1F("mGRefMult", "GRefMult", 500, -3, 3);
   auto * mdTof = new TH1F("#Delta Tof hist", "#Delta Tof", 500, -15, 15);
   auto * mdTofexp = new TH1F("#Delta TofExp hist", "#Delta TofExp", 500, -15, 15);
   auto * mddTof = new TH1F("#Delta #Delta Tof hist", "#Delta #Delta Tof", 500, -10, 10);
   //auto * mZDCEast = new TH1F("mZDCEast", "ZDCEast", 500, 0, 500);
   //auto * mZDCWest = new TH1F("mSDCWest", "ZDCWest", 500, 0, 500);

   //auto * mChargeSum = new TH1F("ChargeSum", "Charge Sum", 500, -3, 3);
   auto * mPt = new TH1F("mPt", "Parent Transverse Momentum", 500, 0, 3);
   //auto * mPhi = new TH1F("mPhi", "Parent Polar Angle (rad)", 500, -3.15, 3.15);
   //auto * mEta = new TH1F("mEta", "Parent Pseudorapidity", 500, -6, 6);
   auto * mMass = new TH1F("mMass", "Parent Mass (GeV)", 500, 0, 2);
   
    // Open the file containing the tree.
    TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    while (myReader.Next()) {

        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e8;
        double me2 = pow(0.00051,2);
        Float_t p1 = pair->d1_mPt;
        Float_t p2 = pair->d2_mPt;
        Float_t p1_2 = pow(p1,2);
        Float_t p2_2 = pow(p2,2);
        Float_t mVertexZVal = pair->mVertexZ;
        UShort_t mGRefMultVal = pair->mGRefMult;  
        Float_t Tof1 = pair->d1_mTof;     
        Float_t Tof2 = pair->d2_mTof;  
        Float_t len1 = pair->d1_mLength;
        Float_t len2 = pair->d2_mLength;    
        Float_t mPtVal = pair->mPt;
        //Float_t mMassVal = pair->mMass;

        //calculate ddTof
        Float_t dTofVal = Tof2 - Tof1;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp2 - texp1;
        Float_t ddTofVal = dTofVal - dTofexpVal;

        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        

        mdTof->Fill( dTofVal );
        mdTofexp->Fill( dTofexpVal );
        mddTof->Fill( ddTofVal );
        mVertexZ->Fill( mVertexZVal);

        //Let's apply some event variable and track cuts
        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4) {
            //now for some track cuts
            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {
                mPt->Fill( lv.Pt() );
                mMass->Fill( lv.M() );
            
            }
        }

        mGRefMult->Fill( mGRefMultVal );




    }
    //fo -> cd();
    
    //makeCanvas();
    //mGRefMult->SetLineColor(kBlack);
    //mGRefMult->Draw();
    //gPad->Print( "plot_mGRefMult.pdf" );

    makeCanvas();
    mPt->SetLineColor(kBlack);
    mPt->GetXaxis()->SetTitle("Pair Transverse Momentum (GeV/c)");
    mPt->GetYaxis()->SetTitle("Counts");
    mPt->Draw();
    gPad->Print( "plot_mPt.pdf" );

    makeCanvas();
    mVertexZ->SetLineColor(kBlack);
    mVertexZ->GetXaxis()->SetTitle("zVertex (cm)");
    mVertexZ->GetYaxis()->SetTitle("Counts");
    mVertexZ->Draw();
    gPad->Print( "plot_mVertexZ.pdf" );

    makeCanvas();
    mMass->SetLineColor(kBlack);
    mMass->GetXaxis()->SetTitle("Parent Mass (GeV/c)");
    mMass->GetYaxis()->SetTitle("Counts");
    mMass->Draw();
    gPad->Print( "plot_mMass.pdf" );

    makeCanvas();
    mdTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mdTof->GetXaxis()->SetTitle("#Delta TOF Distrubition (ns)");
    mdTof->GetYaxis()->SetTitle("Counts");
    mdTof->Draw();
    gPad->Print( "plot_dTof.pdf");

    makeCanvas();
    mdTofexp->SetLineColor(kRed);
    gPad->SetLogy();
    mdTofexp->GetXaxis()->SetTitle("#Delta TOF Distrubition Expected (ns)");
    mdTofexp->GetYaxis()->SetTitle("Counts");
    mdTofexp->Draw();
    gPad->Print( "plot_dTofexp.pdf");

    makeCanvas();
    mddTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mddTof->GetXaxis()->SetTitle("#Delta #Delta TOF Distrubition (ns)");
    mddTof->GetYaxis()->SetTitle("Counts");
    mddTof->Draw();

    //mdTofexp->SetLineColor(kRed);
    //gPad->SetLogy();
    //mdTofexp->Draw();

    gPad->Print( "plot_ddTof.pdf");



    //fo->Write();

}


