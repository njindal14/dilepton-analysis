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
void makeCanvas2() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican6++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.01);

}

void bettereeAnalysis() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   //will instantiate desired histograms below, with numbers describing plots in more detail
   //not including the event variables here, they are used in eventvariables.C

   auto * mPt = new TH1F("mPt", "Parent Transverse Momentum (Gev/c)", 500, 0, 3);
   auto * mPhi = new TH1F("mPhi", "Parent Polar Angle (rad)", 500, -3.15, 3.15);
   auto * mEta = new TH1F("mEta", "Parent Pseudorapidity", 500, -6, 6);
   //auto * mMass = new TH1F("mMass", "Parent Mass (GeV/c^{2})", 500, 0, 4);

   /* Plots to be instiated as follows, referencing ee_note document:
   Currently 5 plots this code makes, all with PID:
   1. Plot the mass distrbution after PID cuts, separated by like sign and unlike sign charge pairs
   2. Plot parent mass without separating (with and without Pt < 100 MeV)
   3. Same as plot 1 but with Pt < 100 MeV
   4. For M > .45 GeV, plot Pt broken up by charge pairs, like or unlike
   5. Chi bands after all PID cuts, recreation of figure 77 in ee_note
    */

   //Plot 1
   auto * masslikecharges = new TH1F("Mass_hist1", "Charge Sum-Separated Parent Mass (GeV) (Plot 1)", 500, 0, 4 );
   auto * massdiffcharges = new TH1F("Mass_hist1", "Charge Sum-Separated Parent Mass (GeV) (Plot 1)", 500, 0, 4);

   //Plot 2
   auto * mMassPt = new TH1F("mMass Pt Cut", "Parent Mass with P_{T} cutoff (GeV) (Plot 2)", 500, 0, 4);
   auto * mMass = new TH1F("mMass", "Parent Mass with P_{T} cutoff (GeV/c^{2}) (Plot 2)", 500, 0, 4); 

   //Plot 3
   auto * Ptmasslikecharges = new TH1F("Mass_hist3", "Charge Sum-Separated Parent Mass (GeV) (Plot 3)", 500, 0, 4 );
   auto * Ptmassdiffcharges = new TH1F("Mass_hist3", "Charge Sum-Separated Parent Mass (GeV) (Plot 3)", 500, 0, 4);

   //Plot 4
   auto * Ptmcutlikesign = new TH1F("Pt_mass_cutoff4", "P_{T} for Parent Mass > 0.45 GeV (Plot 4)", 500, 0, 1);
   auto * Ptmcutdiffsign = new TH1F("Pt_mass_cutoff4", "P_{T} for Parent Mass > 0.45 GeV (Plot 4)", 500, 0, 1);

   //Plot 5
   auto * chiBands2D = new TH2F("Chi Squared Bands", "#chi^{2} Band Comparison (Plot 5)", 200, 0, 200, 200, 0, 160);

   
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

        /* 
        below 15 lines or so calculate necessary things for the PID
        cuts such as Delta Delta TOF, using the math in ee_note document: 
        */
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
        Float_t dTofVal = Tof2 - Tof1;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp2 - texp1;
        Float_t ddTofVal = dTofVal - dTofexpVal;
        
        int chargesumval = pair->mChargeSum;


        
        
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 
        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        /* Now apply necessary PID cuts and fill the histograms:
            Event cuts: Abs val of z vertex < 100, GrefMult <= 4
            then track cuts: abs(DeltaDeltaTof < 0.4), chiee < 10 and 3*chiee < chipipi
        */

 

       if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4) {
            //now for some track cuts
            //For plot 5:
            if(ddTofVal < 0.4 && ddTofVal > -0.4){
                chiBands2D->Fill(chipipi, chiee);   
            }

            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {

                //Conditions for plots 1, 3 and 4:
                if(chargesumval == 0){
                    massdiffcharges->Fill(lv.M());
                    if(lv.Pt() < 0.1){
                        Ptmassdiffcharges->Fill(lv.M());
                    }
                    if(lv.M() > 0.45){
                        Ptmcutdiffsign->Fill(lv.Pt());
                    }
                }
                else{
                    masslikecharges->Fill(lv.M());
                    if(lv.Pt() < 0.1){
                        Ptmasslikecharges->Fill(lv.M());
                    }
                    if(lv.M() > 0.45){
                        Ptmcutlikesign->Fill(lv.Pt());
                    }
                }

                mPt->Fill( lv.Pt() ); 

                //for Plot 2.
                mMass->Fill( lv.M() );
                if(lv.Pt() < 0.1 ){
                    mMassPt->Fill(lv.M());
                }     
            }
        }


     }

     //Make plots

    //Plot 1
    makeCanvas2();

    massdiffcharges->SetLineColor(kBlack);
    gPad->SetLogy();
    massdiffcharges->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    massdiffcharges->GetYaxis()->SetTitle("Counts");
    massdiffcharges->Draw();
    masslikecharges->SetLineColor(kRed);
    masslikecharges->Draw("same");

    TH1F * diff = (TH1F*)massdiffcharges->Clone();

    diff->Add(masslikecharges, -1);
    diff->SetLineColor(kBlue);
    diff->Draw("same");

    auto * leg2 = new TLegend(0.75,0.5,.95,0.7);
    leg2->SetHeader("Legend");
    leg2->AddEntry(masslikecharges,"Like Charges","l");
    leg2->AddEntry(massdiffcharges,"Unlike Charges","l");
    leg2->AddEntry(diff, "Difference", "l");
    leg2->Draw("same");
    gPad->Print( "plot_mMass_by_chargesum.png" );



     //Plot 2
    makeCanvas2();
    mMass->SetLineColor(kBlack);
    mMass->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    mMass->GetYaxis()->SetTitle("Counts");
    mMass->Draw();
    gPad->SetLogy();

    mMassPt->SetLineColor(kRed);
    mMassPt->Draw("same");
    //fit that doesn't work
    mMassPt->Fit("expo");
    //mMasschiTofPt->Fit("fit");

    auto * leg3 = new TLegend(0.75,0.5,.95,0.7);
    leg3->SetHeader("Legend");
    leg3->AddEntry(mMass, "PID cuts","l");
    leg3->AddEntry(mMassPt,"PID cuts and Pt < 100 MeV/c ","l");
    leg3->Draw("same");
    gPad->Print( "plot_mMassPt.png" );


    //Plot 3
    makeCanvas2();
    Ptmassdiffcharges->SetLineColor(kBlack);
    gPad->SetLogy();
    Ptmassdiffcharges->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    Ptmassdiffcharges->GetYaxis()->SetTitle("Counts");
    Ptmassdiffcharges->Draw();
    Ptmasslikecharges->SetLineColor(kRed);
    Ptmasslikecharges->Draw("same");

    TH1F * diff2 = (TH1F*)Ptmassdiffcharges->Clone();

    diff2->Add(Ptmasslikecharges, -1);
    diff2->SetLineColor(kBlue);
    diff2->Draw("same");

    auto * leg4 = new TLegend(0.75,0.5,.95,0.7);
    leg4->SetHeader("Legend");
    leg4->AddEntry(Ptmasslikecharges,"Like Charges, PID + P_{T} < 100 MeV/c","l");
    leg4->AddEntry(Ptmassdiffcharges,"Unlike Charges, PID + P_{T} < 100 MeV/c","l");
    leg4->AddEntry(diff2, "Difference, PID + P_{T} < 100 MeV/c", "l");
    leg4->Draw("same");
    gPad->Print( "plot_mMass_by_chargesum_ptcutoff.png" );

    //Plot 4
    makeCanvas2();
    Ptmcutdiffsign->SetLineColor(kBlack);
    gPad->SetLogy();
    Ptmcutdiffsign->GetXaxis()->SetTitle("Pt (GeV/c)");
    Ptmcutdiffsign->GetYaxis()->SetTitle("Counts");
    Ptmcutdiffsign->Draw();
    Ptmcutlikesign->SetLineColor(kRed);
    Ptmcutlikesign->Draw("same");

    TH1F * diff3 = (TH1F*)Ptmcutdiffsign->Clone();
    diff3->Add(Ptmcutlikesign, -1);
    diff3->SetLineColor(kBlue);
    diff3->Draw("same");

    auto * leg5 = new TLegend(0.75,0.5,.95,0.7);
    leg5->SetHeader("Legend");
    leg5->AddEntry(Ptmcutlikesign,"Like Charges, PID + M > 0.45 GeV/c^{2}","l");
    leg5->AddEntry(Ptmcutdiffsign,"Unlike Charges, PID + M > 0.45 GeV/c^{2}","l");
    leg5->AddEntry(diff3, "Difference, PID + M > 0.45 GeV/c^{2}", "l");
    leg5->Draw("same");
    gPad->Print( "plot_Pt_by_chargesum_mcutoff.png" );

    //Plot 5
    makeCanvas2();
    gPad->SetLogz();
    chiBands2D->SetContour(100);
    chiBands2D->GetXaxis()->SetTitle("#chi_{#pi#pi}^{2}");
    chiBands2D->GetYaxis()->SetTitle("#chi_{ee}^{2}");
    chiBands2D->Draw("colz");




}
