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
    can->SetRightMargin(0.2);

}

double ddToffit(double *x, double *par){
    double fitval;

    fitval = par[0]/(par[2]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + 
            par[3]/(par[5]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[4])/par[5], 2)) +
         + par[6];                   
    return fitval;
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    TLorentzVector lvMinus = lv1 - lv2;
    lv1.Boost(-lvPlus.BoostVector());
    lv2.Boost(-lvPlus.BoostVector());
    double Px = lvPlus.Px();
    double Py = lvPlus.Py();
    double Qx = lvMinus.Px();
    double Qy = lvMinus.Py();
    double absPperp = pow((Px*Px)+(Py*Py), 0.5);
    double absQperp = pow((Qx*Qx)+(Qy*Qy), 0.5);
    double PcrossQ = (Px*Qy) - (Py*Qx);
    double PdotQ = (Px*Qx) + (Py*Qy);
    double cosphi = (Px*Qx + Py*Qy) / (absPperp*absQperp);
    double PairPhi = acos(cosphi);
    if ( PcrossQ > 0 ){
        return PairPhi - 3.141592;
    } else {
        return 3.141592 - PairPhi;
    }
}

void bettereeAnalysis() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   //will instantiate desired histograms below, with numbers describing plots in more detail
   //not including the event variables here, they are used in eventvariables.C

   auto * mPt = new TH1F("mPt", "Parent Transverse Momentum (Gev/c)", 500, 0, 3);
   auto * mPhi = new TH1F("mPhi", "#Delta#phi", 300, -5, 5);
   auto * mEta = new TH1F("mEta", "Parent Pseudorapidity", 500, -6, 6);
   auto * mdTof = new TH1F("#DeltaTOF Hist", "#DeltaTOF", 1000, -15, 15);
   auto * mdTofexp = new TH1F("#DeltaTOFExp Hist", "#DeltaTOFexp", 1000, -15, 15);
   auto * mddTof = new TH1F("#Delta#DeltaTOF Hist", "#Delta#DeltaTOF", 1000, -6, 6);

   auto * ddTofFit = new TF1("fit", ddToffit, -2, 2, 7);
   ddTofFit->SetParameters(100000.0, 0, 0.2, 50000.0, 0, 0.5, 10);
   ddTofFit->SetParNames("A1", "#lambda1", "#sigma1", "A2","#lambda2", "sigma2", "p0");

    //set some parameter ranges
   ddTofFit->SetParLimits(2,0.05, 0.3);
   ddTofFit->SetParLimits(5,0.1,1);

   //smooth it out
   ddTofFit->SetNpx(1000);
   ddTofFit->SetLineWidth(4);

   auto * pathLengths = new TH1F("path lengths", "path lengths", 10000, 0, 500);
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
   auto * mMassPtFit = new TH1F("mMassPtFit", "Fit (Plot 2)", 500, 0,4);

   //Plot 3
   auto * Ptmasslikecharges = new TH1F("Mass_hist3", "Charge Sum-Separated Parent Mass (GeV) (Plot 3)", 500, 0, 4 );
   auto * Ptmassdiffcharges = new TH1F("Mass_hist3", "Charge Sum-Separated Parent Mass (GeV) (Plot 3)", 500, 0, 4);

   //Plot 4
   auto * Ptmcutlikesign = new TH1F("Pt_mass_cutoff4", "P_{T} for Parent Mass > 0.45 GeV (Plot 4)", 500, 0, 1);
   auto * Ptmcutdiffsign = new TH1F("Pt_mass_cutoff4", "P_{T} for Parent Mass > 0.45 GeV (Plot 4)", 500, 0, 1);

   //Plot 5
   auto * chiBands2D = new TH2F("Chi Squared Bands", "#chi^{2} Band Comparison", 200, 0, 200, 200, 0, 160);
   auto * chiBands2DCut = new TH2F("Chi Squared Bands Cut ", "#chi^{2} Band Comparison Cut", 200, 0, 200, 200, 0, 160);


   auto * nSigmaRigidity1 = new TH2F("nsr1", "nsr1", 1000, -2, 2, 1000, -10, 10);
   auto * nSigmaRigidity2 = new TH2F("nsr2", "nsr2", 1000, -2, 2, 1000, -10, 10);

   auto * nSigmaRigidityCut1 = new TH2F("nsrcut1", "nsrcut1", 1000, -2, 2, 1000, -10, 10);
   auto * nSigmaRigidityCut2 = new TH2F("nsrcut2", "nsrcut2", 1000, -2, 2, 1000, -10, 10);



   
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
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;

        /* 
        below 15 lines or so calculate necessary things for the PID
        cuts such as Delta Delta TOF, using the math in ee_note document: 
        */
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
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
        pathLengths->Fill(len1);
        Float_t len2 = pair->d2_mLength; 
        Float_t mPtVal = pair->mPt;
        Float_t dTofVal = Tof1 - Tof2;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp1 - texp2;
        Float_t ddTofVal = dTofVal - dTofexpVal;
        
        int chargesumval = pair->mChargeSum;



        /* Now apply necessary PID cuts and fill the histograms:
            Event cuts: Abs val of z vertex < 100, GrefMult <= 4
            then track cuts: abs(DeltaDeltaTof < 0.4), chiee < 10 and 3*chiee < chipipi
        */
       //double phival = (lv-lvn).Phi();
 

       if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 ) {
            //now for some track cuts
            //For plot 5:
            chiBands2D->Fill(chipipi, chiee); 
            nSigmaRigidity1->Fill(p1, pair->d1_mNSigmaElectron);
            nSigmaRigidity2->Fill(-p2, pair->d2_mNSigmaElectron);

            
            if(ddTofVal < 0.5 && ddTofVal > -0.5){
                // if (true){
                chiBands2DCut->Fill(chipipi, chiee); 
                nSigmaRigidityCut1->Fill(p1, pair->d1_mNSigmaElectron);
                nSigmaRigidityCut2->Fill(-p2, pair->d2_mNSigmaElectron); 

            }

            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi && dca1 < 1 && dca2 < 1) {

                if(lv.M() < 1 && lv.M() > .45) {
                    mPhi->Fill(calc_Phi(lv1,lv2));
                }

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
                    if(lv.M() > 0.5 && lv.M() < 1.5){
                    mMassPtFit->Fill(lv.M());
                }  
                }  
                 
            }
        }
        if(dTofVal != 0 && dTofexpVal !=0 && ddTofVal !=0){
            mdTof->Fill( dTofVal );
            mdTofexp->Fill( dTofexpVal );
            mddTof->Fill( ddTofVal );
        }


     }

     //Make plots

    /*Plot 1
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
    mMassPtFit->SetLineColor(kBlue);
    mMassPtFit->Draw("same");
    mMassPtFit->Fit("expo");
    //mMasschiTofPt->Fit("fit");

    auto * leg3 = new TLegend(0.75,0.5,.95,0.7);
    leg3->SetHeader("Legend");
    leg3->AddEntry(mMass, "PID cuts","l");
    leg3->AddEntry(mMassPt,"PID cuts and Pt < 100 MeV/c ","l");
    leg3->AddEntry(mMassPtFit, "PID cuts, Pt < 100 MeV/c, only section for Fit", "l");
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
    */

    //Plot 5
    makeCanvas2();
    gPad->SetLogz();
    chiBands2D->SetContour(100);
    chiBands2D->GetXaxis()->SetTitle("#chi_{#pi#pi}^{2}");
    chiBands2D->GetYaxis()->SetTitle("#chi_{ee}^{2}");
    chiBands2D->Draw("colz");
    gStyle->SetPalette(1);
    //gPad->SetRightMargin(0.1);
    gPad->Print("plot_chibands.png");

    makeCanvas2();
    gPad->SetLogz();
    chiBands2DCut->SetContour(100);
    chiBands2DCut->GetXaxis()->SetTitle("#chi_{#pi#pi}^{2}");
    chiBands2DCut->GetYaxis()->SetTitle("#chi_{ee}^{2}");
    chiBands2DCut->Draw("colz");
    gStyle->SetPalette(1);
    //gPad->SetRightMargin(0.1);
    gPad->Print("plot_chibandscut.png");
    

    //rigidity plots
    makeCanvas2();
    gPad->SetLogz();
    nSigmaRigidity1->SetContour(100);
    nSigmaRigidity1->GetXaxis()->SetTitle("q*p");
    nSigmaRigidity1->GetYaxis()->SetTitle("n#sigma_{e}1");
    TH2F * sum = (TH2F*)nSigmaRigidity1->Clone();
    sum->Add(nSigmaRigidity2, 1);
    sum->Draw("colz");
    gPad->Print("plot_rigidity.png");


    makeCanvas2();
    gPad->SetLogz();
    nSigmaRigidityCut1->SetContour(100);
    nSigmaRigidityCut1->GetXaxis()->SetTitle("q*p");
    nSigmaRigidityCut1->GetYaxis()->SetTitle("n#sigma_{e}1");
    TH2F * cutsum = (TH2F*)nSigmaRigidityCut1->Clone();
    cutsum->Add(nSigmaRigidityCut2, 1);
    cutsum->Draw("colz");
    gPad->Print("plot_rigiditytofcut.png");

   
    makeCanvas2();
    mdTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mdTof->GetXaxis()->SetTitle("#Delta TOF Distrubition (ns)");
    mdTof->GetYaxis()->SetTitle("Counts");
    mdTof->Draw();
    //gPad->Print( "plot_dTof.png");

    
    //makeCanvas2();
    mdTofexp->SetLineColor(kRed);
    gPad->SetLogy();
    //mdTofexp->GetXaxis()->SetTitle("#Delta TOF Distribution Expected (ns)");
    //mdTofexp->GetYaxis()->SetTitle("Counts");
    mdTofexp->Draw("same");

    auto * legend = new TLegend(0.77,0.6,.97,0.75);
    legend->SetHeader("Legend");
    legend->AddEntry(mdTof,"#DeltaTOF","l");
    legend->AddEntry(mdTofexp,"#DeltaTOFexp","l");
    legend->Draw("same");

    gPad->Print( "plot_dTofPlusExpected.png");

    makeCanvas2();
    mddTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mddTof->GetXaxis()->SetTitle("#Delta #Delta TOF Distrubition (ns)");
    mddTof->GetYaxis()->SetTitle("Counts");
    mddTof->Draw();

    mddTof->Fit("fit", "", "", -2,2);
    gStyle->SetOptFit(1111);
    gPad->Print( "plot_ddTof.png");

    makeCanvas2();
    pathLengths->SetLineColor(kBlack);
    pathLengths->Draw();
    gPad->Print("plot_pathlengths1.png");

    makeCanvas2();
    mPhi->SetLineColor(kBlack);
    mPhi->Draw();
    gPad->Print("plot_phi.png");


    




}
