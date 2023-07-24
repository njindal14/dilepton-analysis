// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH3F.h"

#include "FemtoPairFormat.h"

ClassImp( FemtoPair );

int ican = 0;
void makeCan() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);
}
void rhopAu() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   TFile * fo = new TFile( "data.root", "RECREATE" );

   auto * hulsmass = new TH1F("hulsmass", ";mass (GeV/c^{2})", 500, 0, 3);
   auto * hlspmass = new TH1F("hlspmass", "", 500, 0, 5);
   auto * hlsnmass = new TH1F("hlsnmass", "", 500, 0, 5);

   auto * hulspt = new TH1F("hulspt", ";p_{T} (GeV/c)", 500, 0, 1);
   auto * hlsppt = new TH1F("hlsppt", "", 500, 0, 1);
   auto * hlsnpt = new TH1F("hlsnpt", "", 500, 0, 1);

   auto * hulst = new TH1F("hulst", "p_{T}^{2} (GeV/c)^{2}", 500, 0, 3);
   auto * hlspt = new TH1F("hlspt", "", 500, 0, 5);
   auto * hlsnt = new TH1F("hlsnt", "", 500, 0, 5);

   auto * hpolptmass = new TH3F( "hpolptmass", ";mass;pt;pol", 500, 0, 5, 500, 0, 5, 100, -TMath::Pi(), TMath::Pi() );
   auto * hlspolptmass = new TH3F( "hlspolptmass", ";mass;pt;pol", 500, 0, 5, 500, 0, 5, 100, -TMath::Pi(), TMath::Pi() );

   auto * hc2ptmass = new TH3F( "c2ptmass", ";mass;pt;c2", 500, 0, 5,1000, 0, 1, 100, -1, 1 );
   auto * hlsc2ptmass = new TH3F( "lsc2ptmass", ";mass;pt;c2", 500, 0, 5, 1000, 0, 1, 100, -1, 1 );

   // Open the file containing the tree.
   TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
   
   TTreeReader myReader("PairDst", myFile);

   TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

   
   // TTreeReaderValue<Int_t> qq(myReader, "Pairs.mChargeSum");
   // TTreeReaderValue<Double_t> pt1(myReader, "Pairs.d1_mPt");
   // TTreeReaderValue<Double_t> eta1(myReader, "Pairs.d1_mEta");
   // TTreeReaderValue<Double_t> phi1(myReader, "Pairs.d1_mPhi");
   // TTreeReaderValue<Double_t> dcar1(myReader, "Pairs.d1_mDCA");

   // TTreeReaderValue<Double_t> pt2(myReader, "Pairs.d2_mPt");
   // TTreeReaderValue<Double_t> eta2(myReader, "Pairs.d2_mEta");
   // TTreeReaderValue<Double_t> phi2(myReader, "Pairs.d2_mPhi");
   // TTreeReaderValue<Double_t> dcar2(myReader, "Pairs.d2_mDCA");

   // TTreeReaderValue<Double_t> nSigPi1(myReader, "Pairs.d1_mNSigmaPion");
   // TTreeReaderValue<Double_t> nSigPi2(myReader, "Pairs.d2_mNSigmaPion");
   // // The branch "py" contains floats, too; access those as myPy.
   // // TTreeReaderValue<Double_t> myPy(myReader, "py");
   // // Loop over all entries of the TTree or TChain.
   TLorentzVector lv1, lv2, lv, lvn;
   double h=4;
    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2 ) + pow( pair->d1_mNSigmaPion, 2 );
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
       if ( chipipi < 10 && dca1 < 1 && dca2 < 1){
            
            lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
            lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
            

            lv = lv1 + lv2;
            lvn = lv1 - lv2;
            // if ( gRandom->Rndm() < 0.5 ){
            //     lvn = lv2 - lv1;
            // }

            double pol = lvn.DeltaPhi( lv );
            int qs = pair->mChargeSum;
            TH1F * hmass = hulsmass;
            TH1F * hpt = hulspt;
            TH1F * ht = hulst;
            if ( qs == 0) {
                hmass=hulsmass;
                hpt = hulspt;
                ht = hulst;

                hpolptmass->Fill( lv.M(), lv.Pt(), pol );
                hc2ptmass->Fill( lv.M(), lv.Pt(), cos( h * pol ) );
            }
            else if ( qs == 2){
                hmass = hlspmass;
                hpt = hlsppt;
                ht = hlspt;

                hlspolptmass->Fill( lv.M(), lv.Pt(), pol );
                hlsc2ptmass->Fill( lv.M(), lv.Pt(), cos( h * pol ) );
            }
            else if ( qs == -2) {
                hmass = hlsnmass;
                hpt = hlsnpt;
                ht = hlsnt;

                hlspolptmass->Fill( lv.M(), lv.Pt(), pol );
                hlsc2ptmass->Fill( lv.M(), lv.Pt(), cos( h * pol ) );
            }
            

            hmass->Fill( lv.M() );

            if ( lv.M() > 0.65 && lv.M() < 0.9 ){
                hpt->Fill( lv.Pt() );
                ht->Fill( lv.Pt()*lv.Pt() );
            }
            


       }
   }

   fo->cd();

   makeCan();
   hulsmass->SetLineColor(kBlack);
   hulsmass->Draw();
   hlspmass->SetLineColor(kRed);
   hlsnmass->SetLineColor(kBlue);
   hlsnmass->Draw("same");
   // hlspmass->Draw("same");
   TH1F * hlsmass = (TH1F*)hlspmass->Clone( "hlsmass" );
   hlsmass->Add( hlsnmass );
   hlsmass->SetLineColor(kRed);
   hlsmass->Draw("same");
   TH1F * hsigmass = (TH1F*)hulsmass->Clone( "hsigmass" );
   hsigmass->Add( hlsmass, -1 );
   hsigmass->Draw("same");
   gPad->Print( "plot_m.pdf" );

   makeCan();
   hulspt->SetLineColor(kBlack);
   hulspt->Draw();
   hlsppt->SetLineColor(kRed);
   hlsnpt->SetLineColor(kBlue);
   hlsnpt->Draw("same");
   // hlsppt->Draw("same");
   TH1F * hhlspt = (TH1F*)hlsppt->Clone( "hlspt" );
   hhlspt->Add( hlsnpt );
   hhlspt->SetLineColor(kRed);
   hhlspt->Draw("same");
   TH1F * hsigpt = (TH1F*)hulspt->Clone( "hsigpt" );
   hsigpt->Add( hhlspt, -1 );
   hsigpt->Draw("same");
   gPad->Print( "plot_pt.pdf" );

   makeCan();
   hulst->SetLineColor(kBlack);
   hulst->Draw();
   hlspt->SetLineColor(kRed);
   hlsnt->SetLineColor(kBlue);
   hlsnt->Draw("same");
   // hlspt->Draw("same");
   TH1F * hlst = (TH1F*)hlspt->Clone( "hlst" );
   hlst->Add( hlsnt );
   hlst->SetLineColor(kRed);
   hlst->Draw("same");
   TH1F * hsigt = (TH1F*)hulst->Clone( "hsigt" );
   hsigt->Add( hlst, -1 );
   hsigt->Draw("same");
    gPad->Print( "plot_t.pdf" );


    makeCan();
    hpolptmass->GetYaxis()->SetRangeUser( 0.0, 0.3 );
    hpolptmass->GetXaxis()->SetRangeUser( 0.7, 0.9 );
    hpolptmass->Project3D("z")->Draw("");

    hlspolptmass->GetYaxis()->SetRangeUser( 0.0, 0.3 );
    hlspolptmass->GetXaxis()->SetRangeUser( 0.7, 0.9 );
    hlspolptmass->SetLineColor( TColor::GetColor( "#880e4f" ));
    hlspolptmass->Project3D("z")->Draw("same");


   fo->Write();
  }             
