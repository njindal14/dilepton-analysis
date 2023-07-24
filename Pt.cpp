//includes
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"
#include "TChain.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TH2D.h"



//for plotting
int ican7 = 0;
void makeCanvas4() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican7++ ), "", 900, 600 );
    can->SetTopMargin(0.1);
    can->SetRightMargin(0.05);
    can->SetBottomMargin(0.1);

}

double ddToffit(double *x, double *par){
    double fitval;

    fitval = par[0]/(par[2]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[1])/par[2], 2)) + 
            par[3]/(par[5]* sqrt(2*M_PI)) * exp(-0.5*pow((x[0]-par[4])/par[5], 2)) +
         + par[6];                   
    return fitval;
}

double phiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*(1+par[1]*cos(x[0]) + par[2]*cos(2*x[0]) + par[3]*cos(3*x[0]) + par[4]*cos(4*x[0]));
    return fitval;
}

double chiFit(double *x, double *par){
    double fitval;
    fitval = par[0]*exp(x[0]/par[1]);
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


int main(){

    auto * mPtAu = new TH1F("pT, PID, 0.5 < M_{ee} < 0.8 GeV/c + TOF Cuts, Neutron Selection", "Pair p_{T} (GeV/c)", 50, 0, 0.2);
    auto * mPt = new TH1F("pT, PID, 0.5 < M_{ee} < 0.8 GeV/c + TOF Cuts, Neutron Selection", "Pair p_{T} (GeV/c)", 50, 0, 0.2);
    auto * mPt2 = new TH1F("pT^{2}, PID, 0.5 < M_{ee} < 0.8 GeV/c + TOF Cuts, Neutron Selection", "Pair p_{T}^{2} (GeV/c)", 50, 0, 0.02);


    TFile *QEDPt = new TFile("/Users/Nick/STAR/QEDCode/QED_MB_pt80100_0.40_2.60.root");
    //TH1F * htAu = new TH1F("ht", "Pair p_{T}" , 100, 0, 0.2);
    TH1F * htAu = (TH1F*)QEDPt->Get("ht");   
    
    TFile *QEDPtU = new TFile("/Users/Nick/STAR/QEDCode/QED_MB_pt80100_0.50_0.80Uranium.root");
    //TH1F * htAu = new TH1F("ht", "Pair p_{T}" , 100, 0, 0.2);
    TH1F * htU = (TH1F*)QEDPtU->Get("ht");
    TH1F * htU2 = (TH1F*)QEDPtU->Get("ht2");  




    
    // Open the file containing the tree.

    TFile *myFile = TFile::Open("/Users/Nick/STAR/breit-wheeler/rootFiles/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");


    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run11AuAu.root");
    //TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/slim_pair_dst_Run10AuAu.root");
    TTreeReader myReader2(ch);
    gSystem->Load("FemtoPairFormat_h.so");
    TTreeReaderValue<FemtoPair> pairAu(myReader2, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;


    while (myReader.Next()) {
        //values we will want to use for PID cuts
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        Float_t rapidity = pair->mRapidity;

        /* 
        below 15 lines or so calculate necessary things for the PID
        cuts such as Delta Delta TOF, using the math in ee_note document: 
        */
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        double parentMass = lv.M();
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
        //pathLengths->Fill(len1);
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


        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
       pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {

            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {

                if(lv.M() < 0.8 && lv.M() > .5) {
                    if(mZDCEastVal > 30 && mZDCEastVal < 80 && mZDCWestVal > 30 && mZDCWestVal < 80){
                        mPt->Fill( mPtVal ); 
                        mPt2->Fill(mPtVal*mPtVal);
                    }
                }
            }
       }

    }

    

    while (myReader2.Next()) {
        //values we will want to use for PID cuts
        double chiee = pow( pairAu->d1_mNSigmaElectron, 2 ) + pow( pairAu->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pairAu -> d1_mNSigmaPion, 2) + pow( pairAu -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pairAu->d1_mDCA;
        double dca2 = pairAu->d2_mDCA;
        Float_t rapidity = pairAu->mRapidity;

        /* 
        below 15 lines or so calculate necessary things for the PID
        cuts such as Delta Delta TOF, using the math in ee_note document: 
        */
        
        //Lorentz vectors for each pair track and get lorentz sum and diff
        lv1.SetPtEtaPhiM( pairAu->d1_mPt, pairAu->d1_mEta, pairAu->d1_mPhi, 0.00051 );
        lv2.SetPtEtaPhiM( pairAu->d2_mPt, pairAu->d2_mEta, pairAu->d2_mPhi, 0.00051 ); 

        lv = lv1 + lv2;
        lvn = lv1 - lv2;

        double parentMass = lv.M();
        //get parent total momentum for each track
        Float_t p1 = lv1.P();
        Float_t p2 = lv2.P();
        //square it
        Float_t p1_2 = pow(p1,2);
        Float_t p2_2 = pow(p2,2);
        Float_t mVertexZVal = pairAu->mVertexZ;
        UShort_t mGRefMultVal = pairAu->mGRefMult;  
        Float_t Tof1 = pairAu->d1_mTof;     
        Float_t Tof2 = pairAu->d2_mTof;  
        Float_t len1 = pairAu->d1_mLength;
        //pathLengths->Fill(len1);
        Float_t len2 = pairAu->d2_mLength; 
        Float_t mPtVal = pairAu->mPt;
        Float_t dTofVal = Tof1 - Tof2;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp1 - texp2;
        Float_t ddTofVal = dTofVal - dTofexpVal;
        UShort_t mZDCEastVal = pairAu->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pairAu->mZDCWest;
        
        
        int chargesumval = pairAu->mChargeSum;


        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
       pairAu->d1_mMatchFlag !=0 && pairAu->d2_mMatchFlag!=0) {

            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {

                if(lv.M() < 0.8 && lv.M() > .5) {
                    if(mZDCEastVal > 50 && mZDCEastVal < 350 && mZDCWestVal > 50 && mZDCWestVal < 350){
                        mPtAu->Fill( mPtVal ); 
                    }
                }
            }
        }

    }

    /*gStyle->SetStatBorderSize(0.0);
    gStyle->SetStatStyle(4001);
    gStyle->SetMarkerStyle(1);
    gStyle->SetMarkerSize(5);
    gStyle->SetOptFit(1111);
    gStyle->SetStatX(0.95);
    gStyle->SetStatY(0.5);
    gStyle->SetStatH(0.13);
    gStyle->SetStatW(0.15);*/

    gStyle->SetOptStat(0);

    
    mPtAu->Scale(1/mPtAu->Integral());
    mPt->Scale(1/mPt->Integral());

    htAu->Scale(1/htAu->Integral());
    htU->Scale(1/htU->Integral());
    makeCanvas4();
    mPt->SetLineColor(kBlack);
    mPt->Draw();

    mPtAu->SetLineColor(kRed);
    mPtAu->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    mPtAu->GetYaxis()->SetTitle("Counts");
    mPtAu->Draw("PE;same");
    htAu->SetLineColor(kBlue);
    htAu->Draw("same;hist");
    TLegend * ptlegend = new TLegend(0.75,0.5,.95,0.65);
    ptlegend->AddEntry(mPtAu, "pT Gold", "l");
    ptlegend->AddEntry(htAu,"QED Calculation Gold", "l");
    //ptlegend->Draw("same");
    //gPad->Print("plots/plot_mPtAuQED.png");

    
    
    //htU->Scale(1/htU->Integral());
    //makeCanvas4();
    //mPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    //mPt->GetYaxis()->SetTitle("Counts");
    htU->SetLineColor(kGreen);
    htU->Draw("same;hist");
    //TLegend * ptlegendu = new TLegend(0.79,0.5,.99,.65);
    ptlegend->AddEntry(mPt, "pT Uranium", "l");
    ptlegend->AddEntry(htU, "QED Calculation Uranium", "l");
    //ptlegendu->AddEntry(htU,"QED Calculation", "l");
    ptlegend->Draw("same");
    gPad->Print("plots/plot_mPtBothNeutron1n1nSelectionWithQED.png");

    makeCanvas4();
    TH2D *ptratio = (TH2D*)mPtAu->Clone("ptratio");
    ptratio->SetLineColor(kBlack);
    ptratio->Divide(mPt);
    ptratio->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    ptratio->GetYaxis()->SetTitle("pT Ratio, Au+Au/U+U");
    ptratio->GetYaxis()->SetRangeUser(-2,4);
    ptratio->SetTitle("pT Ratio");
    ptratio->Draw();
    gPad->Print("plots/plot_ptratio.png");

    mPt2->Scale(1/mPt2->Integral());
    htU2->Scale(1/htU2->Integral());
    makeCanvas4();
    mPt2->SetLineColor(kBlack);
    mPt2->Draw("PE");
    htU2->SetLineColor(kGreen);
    htU2->Draw("same;hist");
    gPad->Print("plots/plot_pt2WithQED.png");

    makeCanvas4();
    htU2->Draw("hist");
    gPad->Print("plots/plot_pt2QED.png");


    return 0;


}