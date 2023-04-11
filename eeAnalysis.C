
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
   auto * mZDCEast = new TH1F("mZDCEast", "ZDCEast", 500, 0, 500);
   auto * mZDCWest = new TH1F("mSDCWest", "ZDCWest", 500, 0, 500);

   //auto * mChargeSum = new TH1F("ChargeSum", "Charge Sum", 500, -3, 3);
   auto * mPt = new TH1F("mPt", "Parent Transverse Momentum", 500, 0, 3);
   auto * mPhi = new TH1F("mPhi", "Parent Polar Angle (rad)", 500, -3.15, 3.15);
   auto * mEta = new TH1F("mEta", "Parent Pseudorapidity", 500, -6, 6);
   auto * mMass = new TH1F("mMass", "Parent Mass (GeV)", 500, 0, 4);
   
   
   auto * mMasschi = new TH1F("mMasschi", "Parent Mass (GeV)", 500, 0, 4);
   auto * mMasschiTof = new TH1F("mMasschiTof", "Parent Mass (GeV)", 500, 0, 4);

   
   //new plots

   //for plot 1. indicated below
   auto * masslikecharges = new TH1F("Mass_hist", "Charge Sum-Separated Parent Mass (GeV) (plot 1)", 500, 0, 4 );
   auto * massdiffcharges = new TH1F("Mass_hist", "Charge Sum-Separated Parent Mass (GeV)", 500, 0, 4);

   //for plot 2.1 indicated below
   auto * mMassPt = new TH1F("mMass Pt Cut", "Parent Mass with P_{T} cutoff (GeV) (plot 2.1)", 500, 0, 4);
   auto * mMasschiPt = new TH1F("mMasschi Pt Cut", "Parent Mass (GeV)", 500, 0, 4);
   auto * mMasschiTofPt = new TH1F("mMasschiTof Pt Cut", "Parent Mass (GeV)", 500, 0, 4);
   //auto * func = new TF1("fit", pow(x,[0]), 0, 4);
   //func->SetParameter(0,-4);

   //for plot 2.2 indicated below
   auto * PtcutoffMasslikesign = new TH1F("mass_pt_cutoff", "Parent Mass for P_{T} < 100 MeV (plot 2.2)", 500, 0, 4);
   auto * PtcutoffMassdiffsign = new TH1F("mass_pt_cutoff", "Parent Mass for P_{T} < 100 MeV", 500, 0, 4);
   //auto * massdifference = new TH1F("Mass Difference", "Mass (GeV)", 500, 0, 4);

   //for plot 3 indicated below
   auto * Ptmcutlikesign = new TH1F("Pt_mass_cutoff", "P_{T} for Parent Mass > 0.45 GeV (plot 3)", 500, 0, 1);
   auto * Ptmcutdiffsign = new TH1F("Pt_mass_cutoff", "P_{T} for Parent Mass > 0.45 GeV", 500, 0, 1);

   //recreation of plot 77 from ee analysis note
   auto * chiBands2d = new TH2F("Chi Squared Bands", "#chi^{2} Band Comparison", 200, 0, 200, 200, 0, 160);

   //1. make a plot of unlike-sign pairs and like-sign pairs, also show the difference between the two.
   //2. Make a plot of this plot (shown above) and #1 for a pT < 100 MeV cut, just to show how it changes the phase space distribution.
   //3. Make a plot of the pT distribution for M > 0.45 GeV showing unlike-sign and like-sign (and difference)  
   
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

        //need to calculate different efficiencies and use to find differential cross section
        //how to calculate differential cross section?

        //Figure 77 recreation
        if (ddTofVal < 0.3 && ddTofVal > -0.3) {
            chiBands2d->Fill(chipipi, chiee);  
        }

        mdTof->Fill( dTofVal );
        mdTofexp->Fill( dTofexpVal );
        mddTof->Fill( ddTofVal );
        mVertexZ->Fill( mVertexZVal);
        mMass->Fill( lv.M() );

        int chargesumval = pair->mChargeSum;


        //make plots using like sign and unlike sign pairs for mass, for plot 1. and 2.2
        if(chargesumval == 0){
            massdiffcharges->Fill(lv.M());
            if(lv.Pt() < 0.1){
                PtcutoffMassdiffsign->Fill(lv.M());
            }
        }
        else{
            masslikecharges->Fill(lv.M());
            if(lv.Pt() < 0.1){
                PtcutoffMasslikesign->Fill(lv.M());
            }
        }

        //for plot 3
        if( lv.M() > 0.45 ){
            if(chargesumval == 0){
                Ptmcutdiffsign->Fill(lv.Pt());
            }
            else{
                Ptmcutlikesign->Fill(lv.Pt());
            }
        }

        //make plots for pt < 100 MeV, for plot 2.1
        if(lv.Pt() < .1){
            mMassPt -> Fill(lv.M());
            if(chiee < 5){
                mMasschiPt->Fill( lv.M() );
                if(ddTofVal < 0.5 && ddTofVal > -0.5){
                    mMasschiTofPt->Fill(lv.M());
                }
            }
        }

        if(chiee < 5){
            mMasschi->Fill( lv.M() );
            if(ddTofVal < 0.5 && ddTofVal > -0.5){
                mMasschiTof->Fill(lv.M());              
            }
        }

        //Let's apply some event variable and track cuts
        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4) {
            //now for some track cuts
            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {
                mPt->Fill( lv.Pt() );           
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
    gPad->Print( "plot_mPt.png" );

    /*
    makeCanvas();
    mVertexZ->SetLineColor(kBlack);
    mVertexZ->GetXaxis()->SetTitle("zVertex (cm)");
    mVertexZ->GetYaxis()->SetTitle("Counts");
    mVertexZ->Draw();
    gPad->Print( "plot_mVertexZ.png" );
    */

    makeCanvas();
    //auto * legend = new TLegend(0.2, 0.2, .8, .8);
    mMass->SetLineColor(kBlack);
    mMass->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    mMass->GetYaxis()->SetTitle("Counts");
    mMass->Draw();
    gPad->SetLogy();

    mMasschi->SetLineColor(kRed);
    mMasschi->Draw("same");

    mMasschiTof->SetLineColor(kBlue);
    mMasschiTof->Draw("same");

    auto * leg = new TLegend(0.75,0.5,.95,0.7);
    leg->SetHeader("Legend");
    leg->AddEntry(mMass,"None","l");
    leg->AddEntry(mMasschi,"dE/dx","l");
    leg->AddEntry(mMasschiTof, "dE/dx + TOF", "l");
    leg->Draw("same");
    gPad->Print( "plot_mMass.png" );

    /*
    makeCanvas();
    mdTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mdTof->GetXaxis()->SetTitle("#Delta TOF Distrubition (ns)");
    mdTof->GetYaxis()->SetTitle("Counts");
    mdTof->Draw();
    gPad->Print( "plot_dTof.png");

    
    makeCanvas();
    mdTofexp->SetLineColor(kRed);
    gPad->SetLogy();
    mdTofexp->GetXaxis()->SetTitle("#Delta TOF Distrubition Expected (ns)");
    mdTofexp->GetYaxis()->SetTitle("Counts");
    mdTofexp->Draw();
    gPad->Print( "plot_dTofexp.png");

    makeCanvas();
    mddTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mddTof->GetXaxis()->SetTitle("#Delta #Delta TOF Distrubition (ns)");
    mddTof->GetYaxis()->SetTitle("Counts");
    mddTof->Draw();
    gPad->Print( "plot_ddTof.png");
    */

    //plot 1 as indicated above
    makeCanvas();
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



    //plot 2.1 as indicated above
    makeCanvas();
    mMassPt->SetLineColor(kBlack);
    mMassPt->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    mMassPt->GetYaxis()->SetTitle("Counts");
    mMassPt->Draw();
    gPad->SetLogy();

    mMasschiPt->SetLineColor(kRed);
    mMasschiPt->Draw("same");

    mMasschiTofPt->SetLineColor(kBlue);
    mMasschiTofPt->Draw("same");
    //mMasschiTofPt->Fit("fit");

    auto * leg3 = new TLegend(0.75,0.5,.95,0.7);
    leg3->SetHeader("Legend");
    leg3->AddEntry(mMassPt,"Pt < 100 MeV","l");
    leg3->AddEntry(mMasschiPt,"Pt + dE/dx -- #chi_{ee} < 5","l");
    leg3->AddEntry(mMasschiTofPt, "Pt + dE/dx + TOF", "l");
    leg3->Draw("same");
    gPad->Print( "plot_mMassPt.png" );

    makeCanvas();
    //plot 2.2 as indicated above
    PtcutoffMassdiffsign->SetLineColor(kBlack);
    gPad->SetLogy();
    PtcutoffMassdiffsign->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
    PtcutoffMassdiffsign->GetYaxis()->SetTitle("Counts");
    PtcutoffMassdiffsign->Draw();
    PtcutoffMasslikesign->SetLineColor(kRed);
    PtcutoffMasslikesign->Draw("same");

    TH1F * diff2 = (TH1F*)PtcutoffMassdiffsign->Clone();

    diff2->Add(PtcutoffMasslikesign, -1);
    diff2->SetLineColor(kBlue);
    diff2->Draw("same");

    auto * leg4 = new TLegend(0.75,0.5,.95,0.7);
    leg4->SetHeader("Legend");
    leg4->AddEntry(PtcutoffMasslikesign,"Like Charges","l");
    leg4->AddEntry(PtcutoffMassdiffsign,"Unlike Charges","l");
    leg4->AddEntry(diff2, "Difference", "l");
    leg4->Draw("same");
    gPad->Print( "plot_mMass_by_chargesum_ptcutoff.png" );


    makeCanvas();
    //plot 3 as indicated above
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
    leg5->AddEntry(Ptmcutlikesign,"Like Charges","l");
    leg5->AddEntry(Ptmcutdiffsign,"Unlike Charges","l");
    leg5->AddEntry(diff3, "Difference", "l");
    leg5->Draw("same");
    gPad->Print( "plot_Pt_by_chargesum_mcutoff.png" );

    //Figure 77
    makeCanvas();
    gPad->SetLogz();
    chiBands2d->SetContour(100);
    chiBands2d->GetXaxis()->SetTitle("#chi_{#pi#pi}^{2}");
    chiBands2d->GetYaxis()->SetTitle("#chi_{ee}^{2}");
    chiBands2d->Draw("colz");

    




    //mdTofexp->SetLineColor(kRed);
    //gPad->SetLogy();
    //mdTofexp->Draw();




    //fo->Write();

}


