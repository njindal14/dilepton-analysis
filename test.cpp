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



//for plotting
int ican7 = 0;
void makeCanvas4() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican7++ ), "", 900, 600 );
    can->SetTopMargin(0.04);
    can->SetRightMargin(0.32);

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

    auto * mPt = new TH1F("mPt, PID, M_{ee} + TOF Cuts Au+Au", "Parent Transverse Momentum (GeV/c) Au+Au", 500, 0, 1);
    
    // Open the file containing the tree.
    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run10AuAu.root");
    ch->Add("/Users/Nick/STAR/breit-wheeler/rootFiles/slim_pair_dst_Run11AuAu.root");
    //TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/slim_pair_dst_Run10AuAu.root");
    TTreeReader myReader(ch);
    gSystem->Load("FemtoPairFormat_h.so");
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

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
        
        int chargesumval = pair->mChargeSum;


        if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
       pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {

            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {

                //if(lv.M() < 0.8 && lv.M() > .5) {
                    double phival = calc_Phi(lv1,lv2);
                    mPt->Fill( mPtVal ); 
                //}
            }
       }

    }

    makeCanvas4();
    mPt->SetLineColor(kRed);
    mPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    mPt->GetYaxis()->SetTitle("Counts");
    mPt->Draw();
    gPad->Print("plots/plot_mPtAu.png");

    return 0;


}