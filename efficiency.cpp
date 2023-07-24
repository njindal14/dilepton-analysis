//includes
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TH2F.h"
#include "FemtoPairFormat.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom3.h"
#include <iostream>


//for plotting
int ican8 = 0;
void makeCanvas8() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican8++ ), "", 900, 600 );
    can->SetTopMargin(0.1);
    can->SetBottomMargin(0.1);
    can->SetRightMargin(0.05);
}

double chiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*exp(x[0]/par[1]);
    return fitval;
}

int main() {

    TH1F("h1", "ntuple", 100, -4, 4);
    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;
    gSystem->Load("FemtoPairFormat_h.so");



    auto * massPt = new TH2F("", "Pair Invariant Mass Vs. pT", 300, 0, 1, 300, 0, 2);
    auto * chieephotonconversion = new TH1F("", "#chi_{ee}^{2}", 50, 0, 30);
    auto * Xee = new TH1F("", "#chi_{ee}^{2}", 50, 0, 30);
    auto * chifit = new TF1("chifit", chiFit, 0,5,2);

    auto * STARphiePlusNegativeRapidity = new TH1F("", "STAR #phi e+", 50, -3.14, 3.14);
    auto * STARphiePlusPositiveRapidity = new TH1F("", "STAR #phi e+", 50, -3.14, 3.14);
    auto * STARphiePlusAllRapidity = new TH1F("", "STAR #phi e+", 50, -3.14, 3.14);

    auto * STARphieMinusNegativeRapidity = new TH1F("", "STAR #phi e-", 50, -3.14, 3.14);
    auto * STARphieMinusPositiveRapidity = new TH1F("", "STAR #phi e-", 50, -3.14, 3.14);
    auto * STARphieMinusAllRapidity = new TH1F("", "STAR #phi e-", 50, -3.14, 3.14);

    auto * STARphiPairNegativeRapidity = new TH1F("", "STAR #phi Pair", 50, -3.14, 3.14);
    auto * STARphiPairPositiveRapidity = new TH1F("", "STAR #phi Pair", 50, -3.14, 3.14);
    auto * STARphiPairAllRapidity = new TH1F("", "STAR #phi Pair", 50, -3.14, 3.14);

    auto * phiEtaeMinus = new TH2F("", "Electrons", 50, -1, 1, 50, -3, 3);
    auto * phiEtaePlus = new TH2F("", "Positrons", 50, -1, 1, 50, -3, 3);



    chifit->SetNpx(100);
    chifit->SetLineWidth(4);


    double NPassX2 = 0;
    double NPassTOF = 0;


    while(myReader.Next()){
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

        //lv1 is positron, lv2 is electron

        lv = lv1 + lv2;
        double phi1 = lv1.Phi();
        double phi2 = lv2.Phi();
        double pairPhi = lv.Phi();

        lvn = lv1 - lv2;
        //double parentMass = lv.M();
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
        Float_t len2 = pair->d2_mLength; 
        Float_t mPtVal = pair->mPt;
        Float_t dTofVal = Tof1 - Tof2;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp1 - texp2;
        Float_t ddTofVal = dTofVal - dTofexpVal;  
        int chargesumval = pair->mChargeSum;
       
        //apply cuts for events and tracks (select for BW events)
        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 3 && dca2 < 3) {
            //now for some track cuts
            if( ddTofVal < 0.4 && ddTofVal > -0.4 && ddTofVal !=0) {
                // && chiee < 10 && 3*chiee < chipipi
                massPt->Fill(lv.M(), mPtVal);
                if(lv.M() < 0.02 && mPtVal > 0.3 && mPtVal < 1.2){
                    chieephotonconversion->Fill(chiee);
                }
            }
        }
        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {
            if(ddTofVal < 0.4 && ddTofVal > -0.4 && ddTofVal !=0){
                Xee->Fill(chiee); 
            }
            if(chiee < 10 && 3*chiee < chipipi){
                NPassX2+=1;
                if(ddTofVal < 0.4 && ddTofVal > -0.4 && ddTofVal !=0){
                    NPassTOF+=1;
                }
            }
        }

        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {
            if(chiee < 10 && 3*chiee < chipipi && ddTofVal < 0.4 && ddTofVal > -0.4 && ddTofVal !=0){
                //positron
                STARphiePlusAllRapidity->Fill(phi1);
                STARphieMinusAllRapidity->Fill(phi2);
                STARphiPairAllRapidity->Fill(pairPhi);
                if(lv1.Rapidity() < 0){
                    STARphiePlusNegativeRapidity->Fill(lv1.Phi());
                }
                if(lv1.Rapidity() > 0){
                    STARphiePlusPositiveRapidity->Fill(lv1.Phi());
                }
                //electron
                if(lv2.Rapidity() < 0){
                    STARphieMinusNegativeRapidity->Fill(lv2.Phi());
                    //std::cout << "Filled e minus negatuve rapidity" << lv2.Phi();
                }
                if(lv2.Rapidity() > 0){
                    STARphieMinusPositiveRapidity->Fill(lv2.Phi());
                }
                //pair
                if(lv.Rapidity() < 0){
                    STARphiPairNegativeRapidity->Fill(lv.Phi());
                }
                if(lv.Rapidity() > 0){
                    STARphiPairPositiveRapidity->Fill(lv.Phi());
                }
                phiEtaeMinus->Fill(lv2.Phi(), lv2.Rapidity());
                phiEtaePlus->Fill(lv1.Phi(), lv1.Rapidity());

                
            }
        }
    }



    makeCanvas8();
    massPt->GetXaxis()->SetTitle("Pair Invariant Mass (GeV/c^{2})");
    massPt->GetYaxis()->SetTitle("Pair pT (GeV/c)");
    massPt->SetContour(100);
    massPt->GetXaxis()->SetTitleSize(.05);
    massPt->GetYaxis()->SetTitleSize(.05);
    massPt->Draw("colz");
    gStyle->SetPalette(1);
    gPad->Print("plots/plot_massvpt.png");

    makeCanvas8();
    chieephotonconversion->Scale(1/chieephotonconversion->Integral());
    Xee->Scale(1/Xee->Integral());
    chieephotonconversion->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    chieephotonconversion->GetXaxis()->SetTitleSize(.05);
    chieephotonconversion->GetYaxis()->SetTitleSize(.05);
    chieephotonconversion->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    gPad->SetLogy();
    chieephotonconversion->SetMarkerStyle(20);
    chieephotonconversion->Draw("PE");
    chieephotonconversion->Fit("expo", "", "", 0,5);
    gStyle->SetOptFit(1111);
    chifit->Draw("same");
    //Xee->SetLineColor(kBlue);
    //Xee->Draw("same;PE");
    //TLegend * chiLegend = new TLegend(0.75,0.4,.95,0.55);
    //chiLegend->AddEntry(chieephotonconversion, "Photon Conversion", "l");
    //chiLegend->AddEntry(Xee,"Full Distribution Passing #Delta#Delta TOF", "l");  
    //chiLegend->Draw("same"); 
    gPad->Print("plots/plot_chieePhotonconversionWithFit.png");

    std::cout << "NPassX^2: " << NPassX2;
    std::cout << "NPassTOF: " << NPassTOF;
    std::cout << "TOF Efficiency: " << NPassTOF/NPassX2 * .93 *100 << "%";

    double signalreal = Xee->Integral(Xee->FindFixBin(0), Xee->FindFixBin(10), "");
    double signalideal = chieephotonconversion->Integral(chieephotonconversion->FindFixBin(0), chieephotonconversion->FindFixBin(10), "");
    double efficiencyX2 = signalreal/signalideal;

    int Xeebins = Xee->GetNbinsX();
    double sigPure = 0;
    double sigPureUpper = 0;
    double sigPureLower = 0;
    double sigImpure = 0;
    for(int ix = 1; ix <= Xeebins; ix++){
        if(Xee->GetBinCenter(ix) <= 10){
            sigPure += exp(-1.405-.5591*chieephotonconversion->GetBinCenter(ix));
            sigPureUpper += exp(-1.405+.044 -(.5591-.0236)*chieephotonconversion->GetBinCenter(ix));
            sigPureLower += exp(-1.405-.044 - (.5591+.0236)*chieephotonconversion->GetBinCenter(ix));
            sigImpure += chieephotonconversion->GetBinContent(ix);
        }
    }

    std::cout << "X^2 Efficiency: " << sigPure/sigImpure *100 ;
    std::cout << "X^2 Efficiency lower bound " << sigPureLower/sigImpure *100;
    std::cout << "X^2 Efficiency upper bound " << sigPureUpper/sigImpure *100;


    makeCanvas8();
    STARphieMinusAllRapidity->SetLineColor(kBlack);
    STARphieMinusAllRapidity->GetXaxis()->SetTitle("STAR #phi");
    STARphieMinusAllRapidity->GetYaxis()->SetTitle("Counts");
    STARphieMinusNegativeRapidity->SetLineColor(kRed);
    STARphieMinusPositiveRapidity->SetLineColor(kBlue);
    STARphieMinusAllRapidity->Draw("PE");
    STARphieMinusAllRapidity->GetYaxis()->SetRangeUser(0,400);
    STARphieMinusNegativeRapidity->Draw("same;PE");
    STARphieMinusPositiveRapidity->Draw("same;PE");
    auto * eMinusLegend = new TLegend(0.77,0.5,.97,0.65);
    eMinusLegend->SetHeader("Legend");
    eMinusLegend->AddEntry(STARphieMinusAllRapidity,"e- all rapidity","l");
    eMinusLegend->AddEntry(STARphieMinusPositiveRapidity,"e- #eta > 0","l");
    eMinusLegend->AddEntry(STARphieMinusNegativeRapidity,"e- #eta < 0","l");
    eMinusLegend->Draw("same");
    gPad->Print("plots/plot_STARphieMinus.png");

    makeCanvas8();
    STARphiePlusAllRapidity->SetLineColor(kBlack);
    STARphiePlusAllRapidity->GetXaxis()->SetTitle("STAR #phi");
    STARphiePlusAllRapidity->GetYaxis()->SetTitle("Counts");
    STARphiePlusNegativeRapidity->SetLineColor(kRed);
    STARphiePlusPositiveRapidity->SetLineColor(kBlue);
    STARphiePlusAllRapidity->Draw("PE");
    STARphiePlusAllRapidity->GetYaxis()->SetRangeUser(0,400);
    STARphiePlusNegativeRapidity->Draw("same;PE");
    STARphiePlusPositiveRapidity->Draw("same;PE");
    auto * ePlusLegend = new TLegend(0.77,0.5,.97,0.65);
    ePlusLegend->SetHeader("Legend");
    ePlusLegend->AddEntry(STARphiePlusAllRapidity,"e+ all rapidity","l");
    ePlusLegend->AddEntry(STARphiePlusPositiveRapidity,"e+ #eta > 0","l");
    ePlusLegend->AddEntry(STARphiePlusNegativeRapidity,"e+ #eta < 0","l");
    ePlusLegend->Draw("same");
    gPad->Print("plots/plot_STARphiePlus.png");

    makeCanvas8();
    STARphiPairAllRapidity->SetLineColor(kBlack);
    STARphiPairAllRapidity->GetXaxis()->SetTitle("STAR #phi");
    STARphiPairAllRapidity->GetYaxis()->SetTitle("Counts");
    STARphiPairNegativeRapidity->SetLineColor(kRed);
    STARphiPairPositiveRapidity->SetLineColor(kBlue);
    STARphiPairAllRapidity->Draw("PE");
    STARphiPairAllRapidity->GetYaxis()->SetRangeUser(0,400);
    STARphiPairNegativeRapidity->Draw("same;PE");
    STARphiPairPositiveRapidity->Draw("same;PE");
    auto * pairLegend = new TLegend(0.77,0.5,.97,0.65);
    pairLegend->SetHeader("Legend");
    pairLegend->AddEntry(STARphiPairAllRapidity,"Pair all rapidity","l");
    pairLegend->AddEntry(STARphiPairPositiveRapidity,"Pair #eta > 0","l");
    pairLegend->AddEntry(STARphiPairNegativeRapidity,"Pair #eta < 0","l");
    pairLegend->Draw("same");
    gPad->Print("plots/plot_STARphiPair.png");


    makeCanvas8();
    gPad->SetLogz();
    phiEtaeMinus->SetContour(100);
    phiEtaeMinus->GetXaxis()->SetTitle("#phi");
    phiEtaeMinus->GetYaxis()->SetTitle("#eta");
    gPad->Print("plots/plot_electronPhiEta.png");

    makeCanvas8();
    gPad->SetLogz();
    phiEtaePlus->SetContour(100);
    phiEtaePlus->GetXaxis()->SetTitle("#phi");
    phiEtaePlus->GetYaxis()->SetTitle("#eta");
    gPad->Print("plots/plot_positronPhiEta.png");




}

