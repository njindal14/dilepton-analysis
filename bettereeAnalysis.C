//includes
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"


//for plotting
int ican6 = 0;
void makeCanvas() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican6++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);

}

void eeAnalysis() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   //will instantiate desired histograms below, with numbers describing plots in more detail
   //not including the event variables here, they are used in eventvariables.C

   auto * mPt = new TH1F("mPt", "Parent Transverse Momentum (Gev/c)", 500, 0, 3);
   auto * mPhi = new TH1F("mPhi", "Parent Polar Angle (rad)", 500, -3.15, 3.15);
   auto * mEta = new TH1F("mEta", "Parent Pseudorapidity", 500, -6, 6);
   auto * mMass = new TH1F("mMass", "Parent Mass (GeV/c^{2})", 500, 0, 4);

   /* Plots to be instiated as follows, referencing ee_note document:
   
   
    */





    // Open the file containing the tree.
    TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;


    //loop through events
     while (myReader.Next()) {

        //values we will want to use for PID cuts
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e8;
        double me2 = pow(0.00051,2);

        
        
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 
        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        /* Now apply necessary PID cuts and fill the histograms:


        */


     }

     //Make plots



}
