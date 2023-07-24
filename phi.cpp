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


//for plotting
int ican7 = 0;
void makeCanvas5() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican7++ ), "", 900, 600 );
    can->SetTopMargin(0.1);
    can->SetBottomMargin(0.1);
    can->SetRightMargin(0.05);
}

double phiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*(1+par[1]*cos(x[0]) + par[2]*cos(2*x[0]) + par[3]*cos(3*x[0]) + par[4]*cos(4*x[0]));
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

int main() {

    TH1F("h1", "ntuple", 100, -4, 4);
    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");
    TLorentzVector lv1, lv2, lv, lvn;
    gSystem->Load("FemtoPairFormat_h.so");

    TFile *QEDPhiU = new TFile("/Users/Nick/STAR/QEDCode/QED_MB_deltaphi6080_0.50_0.80.root");
    TH1F * htPhiU = (TH1F*)QEDPhiU->Get("ht");  



    auto * mPhilowest = new TH1F("", "#Delta#phi, 0 < pT < 0.05 GeV/c, 0.5 < M_{ee} < 0.8 GeV/c", 50, -5, 5);
    auto * mPhimid = new TH1F("", "#Delta#phi, 0.025 < pT < 0.075 GeV/c", 50, -5, 5);
    auto * mPhihighest = new TH1F("", "#Delta#phi, 0.05 < pT < 0.15 GeV/c, 0.5 < M_{ee} < 0.8 GeV/c", 50, -5, 5);
    auto * mPhi75to150 = new TH1F("", "#Delta#phi, 0.075 < pT < 0.15 GeV/c, 0.5 < M_{ee} < 0.8 GeV/c", 50, -5, 5);

    auto * mPhiLessthan150 = new TH1F("", "#Delta#phi, pT < 0.15 GeV/c, 0.5 < M_{ee} < 0.8 GeV/c", 50, -5, 5);


    //auto * mPhi150 = new TH1F("mPhi", "#Delta#phi, PID, Mass + TOF cuts, Pt > .15 GeV/c", 100, -5, 5);
    auto * phifit = new TF1("phifit", phiFit, -3.15,3.15,5);
    auto * mPhiZDCHigh = new TH1F("", "#Delta#phi, #sqrt{ZDCEast^{2}+ZDCWest^{2}} > 250, pT < 0.1 GeV/c", 75, -5, 5);
    auto * mPhiZDCLow = new TH1F("", "#Delta#phi, #sqrt{ZDCEast^{2}+ZDCWest^{2}} < 250, pT < 0.1 GeV/c", 75, -5, 5);






    while(myReader.Next()){
        //values we will want to use for PID cuts
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        //UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        //UShort_t mZDCWestVal = pair->mZDCWest;
        //Float_t rapidity = pair->mRapidity;

        /* 
        below 15 lines or so calculate necessary things for the PID
        cuts such as Delta Delta TOF, using the math in ee_note document: 
        */
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
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
        UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pair->mZDCWest;
        
        int chargesumval = pair->mChargeSum;
       


        
        
        //apply cuts for events and tracks (select for BW events)
        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
        pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {
            //now for some track cuts
            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {
                if(lv.M() < 0.8 && lv.M() > .5) {
                    //double phival = calc_Phi(lv1,lv2);
                    double phival;
                    TRandom3 rng(123);
                    double uniform_double = rng.Uniform(0., 1.);
                    if(uniform_double < 0.5){
                        phival = calc_Phi(lv1,lv2);
                    }
                    else{
                        phival = calc_Phi(lv2,lv1);
                    }

                    if(mPtVal < .075 && mPtVal > 0.025){
                        mPhimid->Fill(phival);
                    }
                    if(mPtVal > .05 && mPtVal < .15){
                        mPhihighest->Fill(phival);
                    }
                    if(mPtVal < .05 && mPtVal > 0){
                        mPhilowest->Fill(phival);
                    }
                    if(mPtVal < .1){
                        if(sqrt(pow(mZDCEastVal,2) + pow(mZDCWestVal,2)) < 250){
                            mPhiZDCLow->Fill(phival);
                        }
                        else{
                            mPhiZDCHigh->Fill(phival);
                        }
                    }
                    if(mPtVal < .15 && mPtVal > .075){
                        mPhi75to150->Fill(phival);
                    }
                    if(mPtVal<.15){
                        mPhiLessthan150->Fill(phival);
                    }
                }

            }

        }

    }

    gStyle->SetStatBorderSize(0.0);
    gStyle->SetStatStyle(4001);
    gStyle->SetMarkerStyle(1);
    gStyle->SetMarkerSize(5);
    gStyle->SetOptFit(1111);
    gStyle->SetStatX(0.95);
    gStyle->SetStatY(0.5);
    gStyle->SetStatH(0.13);
    gStyle->SetStatW(0.15);



    makeCanvas5();
    mPhimid->SetLineColor(kBlack);
    mPhimid->Fit("phifit", "", "", -3.15,3.15);
    mPhimid->GetXaxis()->SetTitle("#phi");
    mPhimid->GetYaxis()->SetTitle("Counts");
    mPhimid->Draw("E");
    gPad->Print("plots/plot_phiPtRangemid.png");

    makeCanvas5();
    mPhihighest->SetLineColor(kBlack);
    mPhihighest->Fit("phifit", "", "", -3.15,3.15);
    //gStyle->SetOptFit(1111);
    mPhihighest->GetXaxis()->SetTitle("#phi");
    mPhihighest->GetYaxis()->SetTitle("Counts");
    mPhihighest->Draw("E");
    gPad->Print("plots/plot_phiPtRangehighest.png");

    makeCanvas5();
    mPhilowest->SetLineColor(kBlack);
    mPhilowest->Fit("phifit", "", "", -3.15,3.15);
    //gStyle->SetOptFit(1111);
    mPhilowest->GetXaxis()->SetTitle("#phi");
    mPhilowest->GetYaxis()->SetTitle("Counts");
    mPhilowest->Draw("E");
    gPad->Print("plots/plot_phiPtRangelowest.png");

    makeCanvas5();
    mPhiZDCHigh->SetLineColor(kBlack);
    mPhiZDCHigh->Fit("phifit", "", "", -3.15,3.15);
    //gStyle->SetOptFit(1111);
    mPhiZDCHigh->GetXaxis()->SetTitle("#phi");
    mPhiZDCHigh->GetYaxis()->SetTitle("Counts");
    mPhiZDCHigh->Draw("E");
    gPad->Print("plots/plot_phizdchigh.png");


    makeCanvas5();
    mPhiZDCLow->SetLineColor(kBlack);
    mPhiZDCLow->Fit("phifit", "", "", -3.15,3.15);
    //gStyle->SetOptFit(1111);
    mPhiZDCLow->GetXaxis()->SetTitle("#phi");
    mPhiZDCLow->GetYaxis()->SetTitle("Counts");
    mPhiZDCLow->Draw("E");
    gPad->Print("plots/plot_phizdclow.png");

    makeCanvas5();
    mPhi75to150->SetLineColor(kBlack);
    mPhi75to150->Fit("phifit", "", "", -3.15,3.15);
    mPhi75to150->GetXaxis()->SetTitle("#phi");
    mPhi75to150->GetYaxis()->SetTitle("Counts");
    mPhi75to150->Draw("E");
    gPad->Print("plots/plot_phi75to150.png");


    /*makeCanvas5();
    //mPhiLessthan150->Scale(1/mPhiLessthan150->Integral());
    //htPhiU->Scale(1/htPhiU->Integral());
    mPhiLessthan150->SetLineColor(kBlack);
    mPhiLessthan150->Draw();
    htPhiU->Draw("same;hist");
    TLegend * philegend = new TLegend(0.75,0.5,.95,0.65);
    philegend->AddEntry(mPhiLessthan150, "pT Gold", "l");
    philegend->AddEntry(htPhiU,"QED Calculation", "l");
    philegend->Draw("same");
    gPad->Print("plots/plot_mphiWithQEDCurve.png");*/

    makeCanvas5();
    htPhiU->Scale(1/htPhiU->Integral());
    htPhiU->Fit("phifit", "", "", 0,3.15);
    htPhiU->GetXaxis()->SetTitle("#phi");
    htPhiU->GetYaxis()->SetTitle("Counts");
    htPhiU->SetTitle("#Delta#phi QED Resummation");
    htPhiU->Draw("hist");
    gPad->Print("plots/plot_phiQEDFit.png");

    makeCanvas5();
    mPhiLessthan150->SetLineColor(kBlack);
    mPhiLessthan150->Fit("phifit", "", "", -3.15,3.15);
    mPhiLessthan150->GetXaxis()->SetTitle("#phi");
    mPhiLessthan150->GetYaxis()->SetTitle("Counts");
    mPhiLessthan150->Draw("E");
    gPad->Print("plots/plot_phiPtlessthan150.png");



    return 0;  

}




