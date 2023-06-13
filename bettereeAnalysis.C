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

double* histMoments( TH2F* hist , int n) {
    int nbinsy = hist->GetNbinsY();
    int nbinsx = hist->GetNbinsX();
    double phivals[nbinsx];
    double weights[nbinsx];
    double *cosnphi_moments = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        auto* onedhist = hist->ProjectionX("1dhist", i, i);
	for (int j = 1; j < nbinsx+1; j++) {
            phivals[j-1] = onedhist->GetXaxis()->GetBinCenter(j);
            weights[j-1] = onedhist->GetBinContent(j);
        }
        cosnphi_moments[i-1] = 0;
        for (int k = 0; k < nbinsx; k++){
            if ( onedhist->GetEntries() > 0) {
                cosnphi_moments[i-1] += weights[k] * n * cos(n*phivals[k]) / (onedhist->GetEntries());
            } else { cosnphi_moments[i-1] += 0; }
        }
    }
    
    return cosnphi_moments;
} 

double* moment_error( TH2F* hist, int n ) {
    int nbinsx = hist->GetNbinsX();
    int nbinsy = hist->GetNbinsY();
    double *moment_errors = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        auto* onedhist = hist->ProjectionX("1dhist", i, i);
        auto* momenthist = new TH1F("momenthist", "ncosnphi moments", nbinsx, -1*n, 1*n);
        for (int j = 1; j < nbinsx+1; j++) {
            for (int k = 0; k < onedhist->GetBinContent(j); k++) {
                momenthist->Fill( n * cos( n * (onedhist->GetXaxis()->GetBinCenter(j))));
            }
        }
        moment_errors[i-1] = momenthist->GetMeanError();
    }
    return moment_errors;
}       

double* histYbins( TH2F* hist) {
    int nbinsy = hist->GetNbinsY();
    double *ybinvals = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        ybinvals[i-1] = hist->GetYaxis()->GetBinCenter(i);
    }
    return ybinvals;
}




void bettereeAnalysis() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   //will instantiate desired histograms below, with numbers describing plots in more detail
   //not including the event variables here, they are used in eventvariables.C

   auto * mPt = new TH1F("mPt, PID, M_{ee} + TOF Cuts", "Parent Transverse Momentum (GeV/c)", 40, 0.06, .18);
   auto * mPtAu = new TH1F("mPtAu, PID, M_{ee} + TOF Cuts Au", "Parent Transverse Momentum (GeV/c) Au", 40, 0.06, .18);

   auto * mRapidity = new TH1F("mRapidity, PID + TOF + Pt < 150 MeV + 0.5 < M < 0.8 GeV", "mRapidity, PID + TOF + Pt < 150 MeV + 0.5 < M < 0.8 GeV", 20, -1, 1 );
   auto * mRapidityAu = new TH1F("mRapidity Au, PID + TOF + Pt < 150 MeV + 0.5 < M < 0.8 GeV", "mRapidity, PID + TOF + Pt < 150 MeV + 0.5 < M < 0.8 GeV Au", 20, -1, 1 );


   auto * mPt2 = new TH1F("mPt^{2}, PID, M_{ee} + TOF Cuts, 0.5 < M < 0.8 GeV/c^{2}", "Parent Transverse Momentum Squared (Gev/c)^{2}", 20, 0, 0.0225);
   auto * mPt2Au = new TH1F("mPt^{2}, PID, M_{ee} + TOF Cuts Au, 0.5 < M < 0.8 GeV/c^{2}", "Parent Transverse Momentum Squared (Gev/c)^{2} Au", 20, 0, 0.0225);

   auto * mPhi = new TH1F("mPhi", "#Delta#phi, PID + TOF cuts, Pt < 0.1", 200, -5, 5);
   auto * mPhiAu = new TH1F("mPhiAu", "#Delta#phi, PID + TOF cuts, Pt < 0.1 Au", 200, -5, 5);

   auto * mcosfourphi = new TH1F("mcos4#phi, PID + TOF cuts", "cos(4#phi)", 300, -5, 5);
   auto * mcosthreephi = new TH1F("mcos3#phi, PID + TOF cuts", "cos(3#phi)", 300, -5, 5);
   auto * mcostwophi = new TH1F("mcos2#phi, PID + TOF cuts", "cos(2#phi)", 300, -5, 5);
   auto * mcosphi = new TH1F("mcos#phi, PID + TOF cuts", "cos(#phi)", 300, -5, 5);
   auto * chieeFit = new TF1("chieefit", "[0]", 0, 30);
   chieeFit->SetParameter(0,40);

   auto * PMass = new TH1F("Parent Mass PID + TOF cuts, P_{T} < 0.15 GeV/c", "Parent Mass (GeV/c^{2})", 40, 0, 3);
   auto * PMassAu = new TH1F("Parent Mass PID + TOF cuts Au, P_{T} < 0.15 GeV/c", "Parent Mass (GeV/c^{2})", 40, 0, 3);

   auto * mEta = new TH1F("mEta", "Parent Pseudorapidity", 500, -6, 6);
   auto * mdTof = new TH1F("#DeltaTOF Hist", "#DeltaTOF", 1000, -15, 15);
   auto * mdTofexp = new TH1F("#DeltaTOFExp Hist", "#DeltaTOFexp", 1000, -15, 15);
   auto * mddTof = new TH1F("#Delta#DeltaTOF Hist", "#Delta#DeltaTOF", 1000, -6, 6);
   auto * Xee = new TH1F("#chi_{ee}^{2} Distribution", "#chi_{ee}^{2} Distribution", 100, 0, 30);
   auto * Xee25 = new TH1F("#chi_{ee} Distribution", "#chi_{ee} Distribution, #chi_{#pi#pi} > 25", 200, 0, 15);
   auto * background25 = new TH1F("Background", "#chi_{ee} Distribution background, #chi_{#pi#pi} > 25", 200, 0, 15);
   auto * Xee20 = new TH1F("#chi_{ee} Distribution", "#chi_{ee} Distribution, #chi_{#pi#pi} > 20", 200, 0, 15);
   auto * background20 = new TH1F("Background", "#chi_{ee} Distribution background, #chi_{#pi#pi} > 20", 200, 0, 15);
   auto * Xee15 = new TH1F("#chi_{ee} Distribution", "#chi_{ee} Distribution, #chi_{#pi#pi} > 15", 200, 0, 15);
   auto * background15 = new TH1F("Background", "#chi_{ee} Distribution background, #chi_{#pi#pi} > 15", 200, 0, 15);
   auto * Xee10 = new TH1F("#chi_{ee} Distribution", "#chi_{ee} Distribution, #chi_{#pi#pi} > 10", 200, 0, 15);
   auto * background10 = new TH1F("Background", "#chi_{ee} Distribution background, #chi_{#pi#pi} > 10", 200, 0, 15);
   auto * Xee5 = new TH1F("#chi_{ee} Distribution", "#chi_{ee} Distribution, #chi_{#pi#pi} > 5", 200, 0, 15);
   auto * background5 = new TH1F("Background", "#chi_{ee} Distribution background, #chi_{#pi#pi} > 5", 200, 0, 15);
   auto * Xee1 = new TH1F("#chi_{ee} Distribution", "#chi_{ee} Distribution, #chi_{#pi#pi} > 1", 200, 0, 15);
   auto * background1 = new TH1F("Background", "#chi_{ee} Distribution background, #chi_{#pi#pi} > 1", 200, 0, 15);

   auto * ddTofFit = new TF1("fit", ddToffit, -2, 2, 7);
   auto * phifit = new TF1("phifit", phiFit, -3.15,3.15,5);
   auto * phifitAu = new TF1("phifitAu", phiFit, -3.15,3.15,5);

   auto* chifit = new TF1("chifit", chiFit, 0, 10, 2);
   ddTofFit->SetParameters(100000.0, 0, 0.2, 50000.0, 0, 0.5, 10);
   ddTofFit->SetParNames("A1", "#lambda1", "#sigma1", "A2","#lambda2", "sigma2", "p0");
   phifit->SetParNames("a0", "a1", "a2", "a3", "a4");
   phifit->SetNpx(1000);
   phifit->SetLineWidth(4);
   phifitAu->SetParNames("a0", "a1", "a2", "a3", "a4");
   phifitAu->SetNpx(1000);
   phifitAu->SetLineWidth(4);
   chifit->SetNpx(100);
   chifit->SetLineWidth(4);

    //set some parameter ranges
   ddTofFit->SetParLimits(2,0.05, 0.3);
   ddTofFit->SetParLimits(5,0.1,1);

   //smooth it out
   ddTofFit->SetNpx(1000);
   ddTofFit->SetLineWidth(4);

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


   auto * nSigmaRigidity1 = new TH2F("Rigidity", "Rigidity", 1000, -2, 2, 1000, -10, 10);
   auto * nSigmaRigidity2 = new TH2F("Rigidity2", "Rigidity", 1000, -2, 2, 1000, -10, 10);

   auto * nSigmaRigidityCut1 = new TH2F("Rigidity with #Delta#Delta TOF Cut", "Rigidity with #Delta#Delta TOF Cut", 1000, -2, 2, 1000, -10, 10);
   auto * nSigmaRigidityCut2 = new TH2F("Rigidity2 with #Delta#Delta TOF Cut", "Rigidity with #Delta#Delta TOF Cut", 1000, -2, 2, 1000, -10, 10);

   auto * cos4phivPt = new TH2F("Cos4#phivPt", "cos4#phi distribution vs P_{T}", 400, -2, 2, 15, 0, 0.3);
   auto * cos4phivPtHigherZDC = new TH2F("Cos4#phivPtHigherZDC", "cos4#phi distribution vs P_{T} ADC ZDC Cut", 400, -2, 2, 15, 0, 0.3);
   auto * cos4phivPtLowerZDC = new TH2F("Cos4#phivPtLowerZDC", "cos4#phi distribution vs P_{T} ADC ZDC Cut", 400, -2, 2, 15, 0, 0.3);

   auto * cos3phivPt = new TH2F("Cos3#phivPt", "cos3#phi distribution vs P_{T}", 400, -2, 2, 15, 0, 0.3);
   auto * cos2phivPt = new TH2F("Cos2#phivPt", "cos2#phi distribution vs P_{T}", 400, -2, 2, 15, 0, 0.3);
   auto * cosphivPt = new TH2F("Cos#phivPt", "cos#phi distribution vs P_{T}", 400, -2, 2, 15, 0, 0.3);

   auto * cos4phivPtAu = new TH2F("Cos4#phivPtAu", "cos4#phi distribution vs P_{T} Au", 400, -2, 2, 15, 0, 0.3);
   auto * cos4phivPtHigherZDCAu = new TH2F("Cos4#phivPtHigherZDCAu", "cos4#phi distribution vs P_{T} ADC ZDC Cut Au+Au", 400, -2, 2, 15, 0, 0.3);
   auto * cos4phivPtLowerZDCAu = new TH2F("Cos4#phivPtLowerZDCAu", "cos4#phi distribution vs P_{T} ADC ZDC Cut AU+Au", 400, -2, 2, 15, 0, 0.3);
   auto * cos3phivPtAu = new TH2F("Cos3#phivPtAu", "cos3#phi distribution vs P_{T} Au", 400, -2, 2, 15, 0, 0.3);
   auto * cos2phivPtAu = new TH2F("Cos2#phivPtAu", "cos2#phi distribution vs P_{T} Au", 400, -2, 2, 15, 0, 0.3);
   auto * cosphivPtAu = new TH2F("Cos#phivPtAu", "cos#phi distribution vs P_{T} Au", 400, -2, 2, 15, 0, 0.3);

   auto * cos4phivM = new TH2F("Cos4#phivM", "cos4#phi distribution vs M_{ee}", 100, -2, 2, 15, 0.5, 0.8);
   auto * cos3phivM = new TH2F("Cos3#phivM", "cos3#phi distribution vs M_{ee}", 100, -2, 2, 15, 0.5, 0.8);
   auto * cos2phivM = new TH2F("Cos2#phivM", "cos2#phi distribution vs M_{ee}", 100, -2, 2, 15, 0.5, 0.8);
   auto * cosphivM = new TH2F("Cos#phivM", "cos#phi distribution vs M_{ee}", 100, -2, 2, 15, 0.5, 0.8);

   auto * cos4phivY = new TH2F("Cos4#phivY", "cos4#phi distribution vs Rapidity", 100, -2, 2, 15, -1, 1);
   auto * cos3phivY = new TH2F("Cos3#phivY", "cos3#phi distribution vs Rapidity", 100, -2, 2, 15, -1, 1);
   auto * cos2phivY = new TH2F("Cos2#phivY", "cos2#phi distribution Rapidity", 100, -2, 2, 15, -1, 1);
   auto * cosphivY = new TH2F("Cos#phivY", "cos#phi distribution vs Rapidity", 100, -2, 2, 15, -1, 1);

   auto * phivPt = new TH2F("phivPt", "#phi vs. parent P_{T}; #phi (rad); Parent P_{T} (GeV); Counts", 100, -3.14, 3.14, 100, 0, 0.25);


    // Open the file containing the tree.
    TFile *myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TChain * ch = new TChain("PairDst");
    ch->Add("/Users/Nick/Desktop/Spring2023/slim_pair_dst_Run10AuAu.root");
    ch->Add("/Users/Nick/Desktop/Spring2023/slim_pair_dst_Run11AuAu.root");
    TTreeReader myReader2(ch);
    TTreeReaderValue<FemtoPair> pairAu(myReader2, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;


    //loop through gold-gold data
    while(myReader2.Next()){
        //values we will want to use for PID cuts
        double chiee = pow( pairAu->d1_mNSigmaElectron, 2 ) + pow( pairAu->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pairAu -> d1_mNSigmaPion, 2) + pow( pairAu -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pairAu->d1_mDCA;
        double dca2 = pairAu->d2_mDCA;
        UShort_t mZDCEastVal = pairAu->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pairAu->mZDCWest;
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
        Float_t len2 = pairAu->d2_mLength; 
        Float_t mPtVal = pairAu->mPt;
        Float_t dTofVal = Tof1 - Tof2;
        Float_t texp1 = len1/c * sqrt(1 + me2/p1_2);
        Float_t texp2 = len2/c * sqrt(1 + me2/p2_2);
        Float_t dTofexpVal = texp1 - texp2;
        Float_t ddTofVal = dTofVal - dTofexpVal;
        
        int chargesumval = pairAu->mChargeSum;



       if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
       pairAu->d1_mMatchFlag !=0 && pairAu->d2_mMatchFlag!=0) {
            //now for some track cuts
            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {
                if(mPtVal < 0.15){
                    PMassAu->Fill(lv.M());
                }


                if(lv.M() < 0.8 && lv.M() > .5) {
                    double phival = calc_Phi(lv1,lv2);
                    mPtAu->Fill( mPtVal ); 
                    mPt2Au->Fill( mPtVal * mPtVal);
                    cos4phivPtAu->Fill( 2*cos(4*phival), mPtVal);
                    if(mZDCWestVal > 200 && mZDCEastVal > 200){
                        cos4phivPtHigherZDCAu->Fill(2*cos(4*phival), mPtVal);
                    }
                    else if(mZDCWestVal < 200 && mZDCEastVal < 200){
                        cos4phivPtLowerZDCAu->Fill(2*cos(4*phival), mPtVal);
                    }
                    cos3phivPtAu->Fill( 2*cos(3*phival), mPtVal);
                    cos2phivPtAu->Fill( 2*cos(2*phival), mPtVal);
                    cosphivPtAu->Fill( 2*cos(phival), mPtVal);\

                    if(mPtVal < 0.1){
                        mPhiAu->Fill(phival);
                    }

                    if(mPtVal < .15){
                        mRapidity->Fill(rapidity);
                    }
                }
            }

        }
    }


    //loop through events
     while (myReader.Next()) {

        //values we will want to use for PID cuts
        double chiee = pow( pair->d1_mNSigmaElectron, 2 ) + pow( pair->d2_mNSigmaElectron, 2 );
        double chipipi = pow( pair -> d1_mNSigmaPion, 2) + pow( pair -> d2_mNSigmaPion, 2);
        double c = 3.0e1; //in cm/ns
        double me2 = pow(0.00051,2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        UShort_t mZDCEastVal = pair->mZDCEast;               // ZDC East
        UShort_t mZDCWestVal = pair->mZDCWest;
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
 

       if( mVertexZVal < 100 && mVertexZVal > -100 && mGRefMultVal <= 4 && chargesumval == 0 && dca1 < 1 && dca2 < 1 && 
       pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag!=0) {
            //now for some track cuts
            //For plot 5:
            chiBands2D->Fill(chipipi, chiee); 
            nSigmaRigidity1->Fill(p1, pair->d1_mNSigmaElectron);
            nSigmaRigidity2->Fill(-p2, pair->d2_mNSigmaElectron);



            if(chipipi > 30 &&ddTofVal < 0.4 && ddTofVal > -0.4){
                Xee->Fill(chiee); 
            }
            if(chipipi <= 30 && ddTofVal < 0.4 && ddTofVal > -0.4){

                if(chipipi > 25){
                    if(3*chiee < chipipi ){
                        Xee25->Fill(chiee);
                    }
                    else{
                        background25->Fill(chiee);
                    }
                }
                
                if(chipipi > 20){
                    if(3*chiee < chipipi){
                        Xee20->Fill(chiee);
                    }
                    else{
                        background20->Fill(chiee);
                    }
                }
                
                if(chipipi > 15){
                    if(3*chiee < chipipi){
                        Xee15->Fill(chiee);
                    }
                    else{
                        background15->Fill(chiee);
                    }
                }

                if(chipipi > 10){
                    if(3*chiee < chipipi){
                        Xee10->Fill(chiee);
                    }
                    else{
                        background10->Fill(chiee);
                    }
                }

                if(chipipi > 5){
                    if(3*chiee < chipipi){
                        Xee5->Fill(chiee);
                    }
                    else{
                        background5->Fill(chiee);
                    }
                }
                
                if(chipipi > 1){
                    if(3*chiee < chipipi){
                        Xee1->Fill(chiee);
                    }
                    else{
                        background1->Fill(chiee);
                    }
                }
            }

            
            if(ddTofVal < 0.4 && ddTofVal > -0.4){
                // if (true){
                chiBands2DCut->Fill(chipipi, chiee); 
                nSigmaRigidityCut1->Fill(p1, pair->d1_mNSigmaElectron);
                nSigmaRigidityCut2->Fill(-p2, pair->d2_mNSigmaElectron);
                
            }

            if( ddTofVal < 0.4 && ddTofVal > -0.4 && chiee < 10 && 3*chiee < chipipi) {

                if(mPtVal < 0.15){
                    PMass->Fill(lv.M());
                }

                if(lv.M() < 0.8 && lv.M() > .5) {
                    double phival = calc_Phi(lv1,lv2);
                    mPt->Fill( mPtVal ); 
                    mPt2->Fill( pow(mPtVal , 2));                    
                    phivPt->Fill ( phival, mPtVal);
                    cos4phivPt->Fill( 2*cos(4*phival), mPtVal);
                    cos3phivPt->Fill( 2*cos(3*phival), mPtVal);
                    cos2phivPt->Fill( 2*cos(2*phival), mPtVal);
                    cosphivPt->Fill( 2*cos(phival), mPtVal);
                    if(mZDCWestVal > 200 && mZDCEastVal > 200){
                        cos4phivPtHigherZDC->Fill(2*cos(4*phival), mPtVal);
                    }
                    else if(mZDCWestVal < 200 && mZDCEastVal < 200){
                        cos4phivPtLowerZDC->Fill(2*cos(4*phival), mPtVal);
                    }
                    cos4phivM->Fill( 2*cos(4*phival), parentMass);
                    cos3phivM->Fill( 2*cos(3*phival), parentMass);
                    cos2phivM->Fill( 2*cos(2*phival), parentMass);
                    cosphivM->Fill( 2*cos(phival), parentMass);

                    cos4phivY->Fill( 2*cos(4*phival), rapidity);
                    cos3phivY->Fill( 2*cos(3*phival), rapidity);
                    cos2phivY->Fill( 2*cos(2*phival), rapidity);
                    cosphivY->Fill( 2*cos(phival), rapidity);
                    

                    if(mPtVal < 0.1){
                        mPhi->Fill(phival);
                        mcosfourphi->Fill(cos(4*phival));
                    }
                    if(mPtVal < 0.15){
                        mRapidityAu->Fill(rapidity);
                    }
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
    /*makeCanvas2();
    gPad->SetLogz();
    chiBands2D->SetContour(100);
    chiBands2D->GetXaxis()->SetTitle("#chi_{#pi#pi}^{2}");
    chiBands2D->GetYaxis()->SetTitle("#chi_{ee}^{2}");
    chiBands2D->Draw("colz");
    gStyle->SetPalette(1);
    //gPad->SetRightMargin(0.1);
    gPad->Print("plots/plot_chibands.png");

    makeCanvas2();
    gPad->SetLogz();
    chiBands2DCut->SetContour(100);
    chiBands2DCut->GetXaxis()->SetTitle("#chi_{#pi#pi}^{2}");
    chiBands2DCut->GetYaxis()->SetTitle("#chi_{ee}^{2}");
    chiBands2DCut->Draw("colz");
    gStyle->SetPalette(1);
    //gPad->SetRightMargin(0.1);
    gPad->Print("plots/plot_chibandscut.png");
    

    //rigidity plots
    makeCanvas2();
    gPad->SetLogz();
    nSigmaRigidity1->SetContour(100);
    nSigmaRigidity1->GetXaxis()->SetTitle("q*p");
    nSigmaRigidity1->GetYaxis()->SetTitle("n#sigma_{e}1");
    TH2F * sum = (TH2F*)nSigmaRigidity1->Clone();
    sum->Add(nSigmaRigidity2, 1);
    sum->Draw("colz");
    gPad->Print("plots/plot_rigidity.png");


    makeCanvas2();
    gPad->SetLogz();
    nSigmaRigidityCut1->SetContour(100);
    nSigmaRigidityCut1->GetXaxis()->SetTitle("q*p (GeV/c)");
    nSigmaRigidityCut1->GetYaxis()->SetTitle("n#sigma_{e1}");
    TH2F * cutsum = (TH2F*)nSigmaRigidityCut1->Clone();
    cutsum->Add(nSigmaRigidityCut2, 1);
    cutsum->Draw("colz");
    gPad->Print("plots/plot_rigiditytofcut.png"); 

   
    makeCanvas2();
    mdTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mdTof->GetXaxis()->SetTitle("#Delta TOF Distrubition (ns)");
    mdTof->GetYaxis()->SetTitle("Counts");
    mdTof->Draw();
    gPad->Print( "plots/plot_dTof.png");

    
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

    gPad->Print( "plots/plot_dTofPlusExpected.png");

    makeCanvas2();
    mddTof->SetLineColor(kBlack);
    gPad->SetLogy();
    mddTof->GetXaxis()->SetTitle("#Delta #Delta TOF Distrubition (ns)");
    mddTof->GetYaxis()->SetTitle("Counts");
    mddTof->Draw();

    mddTof->Fit("fit", "", "", -2,2);
    gStyle->SetOptFit(1111);
    gPad->Print( "plots/plot_ddTof.png"); */

    makeCanvas2();
    mPhi->SetLineColor(kBlack);
    mPhi->Fit("phifit", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111);
    mPhi->GetXaxis()->SetTitle("#phi");
    mPhi->GetYaxis()->SetTitle("Counts");
    mPhi->Draw();
    gPad->Print("plots/plot_phi.png");

    makeCanvas2();
    mPhiAu->SetLineColor(kBlack);
    mPhiAu->Fit("phifitAu", "", "", -3.15,3.15);
    gStyle->SetOptFit(1111);
    mPhiAu->GetXaxis()->SetTitle("#phi Au");
    mPhiAu->GetYaxis()->SetTitle("Counts");
    mPhiAu->Draw();
    gPad->Print("plots/plot_phiAu.png");

    makeCanvas2();
    mcosfourphi->SetLineColor(kBlack);
    mcosfourphi->GetXaxis()->SetTitle("cos(4#phi)");
    mcosfourphi->GetYaxis()->SetTitle("Counts");
    mcosfourphi->Draw();
    gPad->Print("plots/plot_cos(4phi).png");


    makeCanvas2();
    cos4phivPt->GetXaxis()->SetTitle("cos(4#phi)");
    cos4phivPt->GetYaxis()->SetTitle("P_{T} (GeV/c)");
    cos4phivPt->Draw("colz");
    gPad->Print( "plots/plot_cos4phivPt.png" );       

    makeCanvas2();
    phivPt->Draw("colz");
    phivPt->GetXaxis()->SetTitle("#phi");
    phivPt->GetYaxis()->SetTitle("P_{T}");
    gPad->Print("plots/plot_phivPt.png");

    auto *m2Ptcos4phimoments = cos4phivPt->ProfileY("m2Ptcos4phimoments",1, -1);
    auto *m2Ptcos3phimoments = cos3phivPt->ProfileY("m2Ptcos3phimoments", 1, -1);
    auto *m2Ptcos2phimoments = cos2phivPt->ProfileY("m2Ptcos2phimoments",1, -1);
    auto *m2Ptcosphimoments = cosphivPt->ProfileY("m2Ptcosphimoments",1, -1);


    auto *m2Ptcos4phimomentsAu = cos4phivPtAu->ProfileY("m2Ptcos4phimomentsAu", 1, -1);
    auto *m2Ptcos3phimomentsAu = cos3phivPtAu->ProfileY("m2Ptcos3phimomentsAu",1, -1);
    auto *m2Ptcos2phimomentsAu = cos2phivPtAu->ProfileY("m2Ptcos2phimomentsAu", 1, -1);
    auto *m2PtcosphimomentsAu = cosphivPtAu->ProfileY("m2PtcosphimomentsAu",1, -1);

    auto *m2Mcos4phimoments = cos4phivM->ProfileY("m2Mcos4phimoments", 1, -1);
    auto *m2Mcos3phimoments = cos3phivM->ProfileY("m2Mcos3phimoments", 1, -1);
    auto *m2Mcos2phimoments = cos2phivM->ProfileY("m2Mcos2phimoments", 1, -1);
    auto *m2Mcosphimoments = cosphivM->ProfileY("m2Mcosphimoments", 1, -1);

    auto *m2Ycos4phimoments = cos4phivY->ProfileY("m2Ycos4phimoments", 1, -1);
    auto *m2Ycos3phimoments = cos3phivY->ProfileY("m2Ycos3phimoments", 1, -1);
    auto *m2Ycos2phimoments = cos2phivY->ProfileY("m2Ycos2phimoments", 1, -1);
    auto *m2Ycosphimoments = cosphivY->ProfileY("m2Ycosphimoments", 1, -1);




    makeCanvas2();
    m2Mcos4phimoments->SetTitle("scaled cos(4#phi) moments vs. M_{ee}; M_{ee} (GeV/c^{2}); 2<cos(4#phi)>");
    //m2Mcos4phimoments->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    //m2Mcos4phimoments->GetYaxis()->SetTitle("cos(4#phi) moments scaled");
    m2Mcos4phimoments->SetLineColor(kBlack);
    m2Mcos4phimoments->Draw();
    gPad->Print( "plots/plot_m2Mcos4phimoments.png");


    makeCanvas2();
    m2Mcos3phimoments->SetTitle("scaled cos(3#phi) moments vs. M_{ee}; M_{ee} (GeV/c^{2}); 2<cos(3#phi)>");
    //m2Mcos3phimoments->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    //m2Mcos3phimoments->GetYaxis()->SetTitle("cos(3#phi) moments scaled");
    m2Mcos3phimoments->SetLineColor(kRed);
    m2Mcos3phimoments->Draw();
    gPad->Print( "plots/plot_m2Mcos3phimoments.png");


    makeCanvas2();
    m2Mcos2phimoments->SetTitle("scaled cos(2#phi) moments vs. M_{ee}; M_{ee} (GeV/c^{2}); 2<cos(2#phi)>");
    //m2Mcos2phimoments->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    //m2Mcos2phimoments->GetYaxis()->SetTitle("cos(2#phi) moments scaled");
    m2Mcos2phimoments->SetLineColor(kBlack);
    m2Mcos2phimoments->Draw();
    gPad->Print( "plots/plot_m2Mcos2phimoments.png");

    makeCanvas2();
    m2Mcosphimoments->SetTitle("scaled cos(#phi) moments vs. M_{ee}; M_{ee} (GeV/c^{2}); 2<cos(#phi)>");
    //m2Mcosphimoments->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    //m2Mcosphimoments->GetYaxis()->SetTitle("cos(#phi) moments scaled");
    m2Mcosphimoments->SetLineColor(kBlack);
    m2Mcosphimoments->Draw();
    gPad->Print( "plots/plot_m2Mcosphimoments.png");



    makeCanvas2();
    m2Ycos4phimoments->SetTitle("Cos(4#phi) Moments vs. Rapidity; Rapidity; 2<cos(4#phi)>");
    m2Ycos4phimoments->SetLineColor(kBlack);
    m2Ycos4phimoments->Draw();
    gPad->Print( "plots/plot_m2Ycos4phimoments.png");


    makeCanvas2();
    m2Ycos3phimoments->SetTitle("Cos(3#phi) Moments vs. Rapidity; Rapidity; 2<cos(3#phi)>");
    m2Ycos3phimoments->SetLineColor(kRed);
    m2Ycos3phimoments->Draw();
    gPad->Print( "plots/plot_m2Ycos3phimoments.png");


    makeCanvas2();
    m2Ycos2phimoments->SetTitle("Cos(2#phi) moments vs. Rapidity; Rapidity; 2<cos(2#phi)>");
    m2Ycos2phimoments->SetLineColor(kBlack);
    m2Ycos2phimoments->Draw();
    gPad->Print( "plots/plot_m2Ycos2phimoments.png");

    makeCanvas2();
    m2Ycosphimoments->SetTitle("Cos(#phi) moments vs. Rapidity; Rapidity; 2<cos(#phi)>");
    m2Ycosphimoments->SetLineColor(kBlack);
    m2Ycosphimoments->Draw();
    gPad->Print( "plots/plot_m2Ycosphimoments.png");

    makeCanvas2();
    Xee->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee->SetMarkerStyle(20);
    Xee->Draw("PE");
    chieeFit->SetLineWidth(4);
    chieeFit->SetLineColor(kBlue);
    Xee->Fit("chieefit", "", "", 12, 30);
    Xee->Fit("expo", "R+", "", 0., 12.);
    gStyle->SetOptFit(1111);
    chieeFit->Draw("same");
    double background0 = chieeFit->GetParameter(0);
    double binsx0 = Xee->GetNbinsX();
    double totbackground = 0;
    double totsignal = 0;
    for(int ix =1; ix <= Xee->FindFixBin(10); ix++){
        totsignal += Xee->GetBinContent(ix);
        totbackground += background0;
    }
    double purity = totsignal/ (totsignal+totbackground);
    cout << "Purity for chipipi > 30: " << purity*100 << "%\n";
    gPad->Print("plots/chi2eePlot.png");

    makeCanvas2();
    Xee25->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee25->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee25->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee25->Draw();
    background25->SetLineColor(kRed);
    background25->Draw("same");
    auto * legend25 = new TLegend(0.77,0.5,.97,0.65);
    legend25->SetHeader("Legend");
    legend25->AddEntry(Xee25,"Signal","l");
    legend25->AddEntry(background25,"Background","l");
    legend25->Draw("same");
    TLine * l = new TLine(8.33,0,8.33,60);
    l->SetLineColor(kGreen);
    l->Draw("same");
    int xbins25 = Xee25->GetNbinsX();
    double sig25 = 0;
    double bground25 = 0;
    for(int ix = 1; ix <= xbins25; ix++){
        if(Xee25->GetBinCenter(ix) <= 8.33){
            sig25 += Xee25->GetBinContent(ix);
            bground25 += background25->GetBinContent(ix);
        }
    }
    cout << "Signal 25: " << sig25 << "\n";
    cout << "Background 25: " << bground25 << "\n";
    double purity25 = sig25/(sig25 + bground25);
    cout << "Purity25 : " << purity25 << "\n";
    double sigint25 = Xee25->Integral(Xee25->FindFixBin(0), Xee25->FindFixBin(8.33), "");
    double bkg25 = background25->Integral(background25->FindFixBin(0), background25->FindFixBin(8.33), "");
    double purityIntegral25 = sigint25/(sigint25+bkg25);
    cout << "Purity 25 Integral Method: " << purityIntegral25;
    gPad->Print("plots/plot_chi25.png");


    makeCanvas2();
    Xee20->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee20->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee20->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee20->Draw();
    background20->SetLineColor(kRed);
    background20->Draw("same");
    auto * legend20 = new TLegend(0.77,0.5,.97,0.65);
    legend20->SetHeader("Legend");
    legend20->AddEntry(Xee20,"Signal","l");
    legend20->AddEntry(background20,"Background","l");
    legend20->Draw("same");
    TLine * l2 = new TLine(6.66,0,6.66,60);
    l2->SetLineColor(kGreen);
    l2->Draw("same");
    
    int xbins20 = Xee20->GetNbinsX();
    double sig20 = 0;
    double bground20 = 0;
    for(int ix = 1; ix <= xbins20; ix++){
        if(Xee20->GetBinCenter(ix) <= 6.66){
            sig20 += Xee20->GetBinContent(ix);
            bground20 += background20->GetBinContent(ix);
        }
    }
    cout << "Signal 20: " << sig20 << "\n";
    cout << "Background 20: " << bground20 << "\n";
    double purity20 = sig20/(sig20 + bground20);
    cout << "Purity20 : " << purity20 << "\n";
    double sigint20 = Xee20->Integral(Xee20->FindFixBin(0), Xee20->FindFixBin(6.66), "");
    double bkg20 = background20->Integral(background20->FindFixBin(0), background20->FindFixBin(6.66), "");
    double purityIntegral20 = sigint20/(sigint20+bkg20);
    cout << "Purity 20 Integral Method: " << purityIntegral20;
    gPad->Print("plots/plot_chi20.png");


    makeCanvas2();
    Xee15->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee15->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee15->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee15->Draw();
    background15->SetLineColor(kRed);
    background15->Draw("same");
    auto * legend15 = new TLegend(0.77,0.5,.97,0.65);
    legend15->SetHeader("Legend");
    legend15->AddEntry(Xee15,"Signal","l");
    legend15->AddEntry(background15,"Background","l");
    legend15->Draw("same");
    TLine * l3 = new TLine(5,0,5,60);
    l3->SetLineColor(kGreen);
    l3->Draw("same");
    int xbins15 = Xee15->GetNbinsX();
    double sig15 = 0;
    double bground15 = 0;
    for(int ix = 1; ix <= xbins15; ix++){
        if(Xee15->GetBinCenter(ix) <= 5){
            sig15 += Xee15->GetBinContent(ix);
            bground15 += background15->GetBinContent(ix);
        }
    }
    cout << "Signal 15: " << sig15 << "\n";
    cout << "Background 15: " << bground15 << "\n";
    double purity15 = sig15/(sig15 + bground15);
    cout << "Purity15 : " << purity15 << "\n";
    double sigint15 = Xee15->Integral(Xee15->FindFixBin(0), Xee15->FindFixBin(5), "");
    double bkg15 = background15->Integral(background15->FindFixBin(0), background15->FindFixBin(5), "");
    double purityIntegral15 = sigint15/(sigint15+bkg15);
    cout << "Purity 15 Integral Method: " << purityIntegral15;
    gPad->Print("plots/plot_chi15.png");


    makeCanvas2();
    Xee10->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee10->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee10->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee10->Draw();
    background10->SetLineColor(kRed);
    background10->Draw("same");
    auto * legend10 = new TLegend(0.77,0.5,.97,0.65);
    legend10->SetHeader("Legend");
    legend10->AddEntry(Xee10,"Signal","l");
    legend10->AddEntry(background10,"Background","l");
    legend10->Draw("same");
    TLine * l4 = new TLine(3.33,0,3.33,60);
    l4->SetLineColor(kGreen);
    l4->Draw("same");
    int xbins10 = Xee10->GetNbinsX();
    double sig10 = 0;
    double bground10 = 0;
    for(int ix = 1; ix <= xbins10; ix++){
        if(Xee10->GetBinCenter(ix) <= 3.33){
            sig10 += Xee10->GetBinContent(ix);
            bground10 += background10->GetBinContent(ix);
        }
    }
    cout << "Signal 10: " << sig10 << "\n";
    cout << "Background 10: " << bground10 << "\n";
    double purity10 = sig10/(sig10 + bground10);
    cout << "Purity10 : " << purity10 << "\n";
    double sigint10 = Xee10->Integral(Xee10->FindFixBin(0), Xee10->FindFixBin(3.33), "");
    double bkg10 = background10->Integral(background10->FindFixBin(0), background10->FindFixBin(3.33), "");
    double purityIntegral10 = sigint10/(sigint10+bkg10);
    cout << "Purity 10 Integral Method: " << purityIntegral10;
    gPad->Print("plots/plot_chi10.png");



    makeCanvas2();
    Xee5->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee5->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee5->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee5->Draw();
    background5->SetLineColor(kRed);
    background5->Draw("same");
    auto * legend5 = new TLegend(0.77,0.5,.97,0.65);
    legend5->SetHeader("Legend");
    legend5->AddEntry(Xee5,"Signal","l");
    legend5->AddEntry(background5,"Background","l");
    legend5->Draw("same");
    TLine * l5 = new TLine(1.66,0,1.66,60);
    l5->SetLineColor(kGreen);
    l5->Draw("same");
    int xbins5 = Xee5->GetNbinsX();
    double sig5 = 0;
    double bground5 = 0;
    for(int ix = 1; ix <= xbins5; ix++){
        if(Xee5->GetBinCenter(ix) <= 1.66){
            sig5 += Xee5->GetBinContent(ix);
            bground5 += background5->GetBinContent(ix);
        }
    }
    cout << "Signal 5: " << sig5 << "\n";
    cout << "Background 5: " << bground5 << "\n";
    double purity5 = sig5/(sig5 + bground5);
    cout << "Purity5 : " << purity5 << "\n";
    double sigint5 = Xee5->Integral(Xee5->FindFixBin(0), Xee5->FindFixBin(1.66), "");
    double bkg5 = background5->Integral(background5->FindFixBin(0), background5->FindFixBin(1.66), "");
    double purityIntegral5 = sigint5/(sigint5+bkg5);
    cout << "Purity 5 Integral Method: " << purityIntegral5;
    gPad->Print("plots/plot_chi5.png");


    makeCanvas2();
    Xee1->SetLineColor(kBlack);
    gPad->SetLogy();
    Xee1->GetXaxis()->SetTitle("#chi_{ee}^{2}");
    Xee1->GetYaxis()->SetTitle("dN/d(#chi_{ee}^{2})");
    Xee1->Draw();
    background1->SetLineColor(kRed);
    background1->Draw("same");
    auto * legend1 = new TLegend(0.77,0.5,.97,0.65);
    legend1->SetHeader("Legend");
    legend1->AddEntry(Xee1,"Signal","l");
    legend1->AddEntry(background1,"Background","l");
    legend1->Draw("same");
    TLine * l6 = new TLine(.33,0,.33,60);
    l6->SetLineColor(kGreen);
    l6->Draw("same");
    int xbins1 = Xee1->GetNbinsX();
    double sig1 = 0;
    double bground1 = 0;
    for(int ix = 1; ix <= xbins1; ix++){
        if(Xee1->GetBinCenter(ix) <= .33){
            sig1 += Xee1->GetBinContent(ix);
            bground1 += background1->GetBinContent(ix);
        }
    }
    cout << "Signal 1: " << sig1 << "\n";
    cout << "Background 1: " << bground1 << "\n";
    double purity1 = sig1/(sig1 + bground1);
    cout << "Purity1 : " << purity1 << "\n";
    double sigint1 = Xee1->Integral(Xee1->FindFixBin(0), Xee1->FindFixBin(.33), "");
    double bkg1 = background1->Integral(background1->FindFixBin(0), background1->FindFixBin(.33), "");
    double purityIntegral1 = sigint1/(sigint1+bkg1);
    cout << "Purity 1 Integral Method: " << purityIntegral1;

    gPad->Print("plots/plot_chi1.png");

    double overallPurity = 0;
    overallPurity += purity1 + purity10 + purity15 + purity5 + purity20 + purity25;
    overallPurity /= 6;
    cout << "AVERAGE SLICE PURITY: " << overallPurity*100 << "%\n";

    double overallPurityIntegral = 0;
    overallPurityIntegral+= purityIntegral1 + purityIntegral5 + purityIntegral10 + purityIntegral15 + purityIntegral20 + purityIntegral25;
    overallPurityIntegral/=6;
    cout << "AVERAGE SLICE PURITY INTEGRAL METHOD: " << overallPurity*100 << "%\n";


    makeCanvas2();
    mPtAu->SetLineColor(kBlack);
    mPtAu->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    mPtAu->GetYaxis()->SetTitle("1/P_{T} * dN/d(P_{T})");
    for(int ix = 1; ix < mPtAu->GetNbinsX(); ix++){
        mPtAu->SetBinContent(ix, 1/mPtAu->GetBinCenter(ix) * mPtAu->GetBinContent(ix));
    }
    //mPtAu->Scale(1/mPtAu);
    mPtAu->Scale(1/mPtAu->GetEntries());
    mPtAu->Draw();
    gPad->Print("plots/plot_mPtAuHighPtRange.png");


    makeCanvas2();
    mPt->SetLineColor(kRed);
    mPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    mPt->GetYaxis()->SetTitle("1/P_{T} * dN/d(P_{T})");
    for(int ix = 1; ix < mPt->GetNbinsX(); ix++){
        mPt->SetBinContent(ix, 1/mPt->GetBinCenter(ix) * mPt->GetBinContent(ix));
    }
    //mPt->Scale(1/mPt);
    mPt->Scale(1/mPt->GetEntries());
    mPt->Draw();
    gPad->Print("plots/plot_mPtHighPtRange.png");

    makeCanvas2();
    TH1F *mPtRatio = (TH1F*)mPt->Clone("mPtRatio");
    mPtRatio->SetLineColor(kBlack);
    mPtRatio->Divide(mPtAu);
    mPtRatio->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    mPtRatio->GetYaxis()->SetTitle("P_{T} Ratio, U+U/Au+Au");
    mPtRatio->GetXaxis()->SetRangeUser(0,0.15);
    mPtRatio->GetYaxis()->SetRangeUser(0,2);
    mPtRatio->SetTitle("P_{T} Ratio");
    mPtRatio->Fit("pol1", "", "", 0.01, 0.1);
    mPtRatio->Draw();
    gPad->Print("plots/plot_mPtRatioUAuHighPtRange.png");

    makeCanvas2();
    PMass->SetLineColor(kBlack);
    PMassAu->SetLineColor(kRed);
    PMass->GetXaxis()->SetTitle("Pair Invariant Mass (GeV/c^{2})");
    PMass->GetYaxis()->SetTitle("Counts");
    PMass->SetTitle("Pair Invariant Mass (GeV/c^{2})");
    PMass->Scale(1/PMass->GetEntries());
    PMassAu->Scale(1/PMassAu->GetEntries());
    PMass->Draw();
    PMassAu->Draw("same");
    auto * legendMass = new TLegend(0.79,0.5,.99,0.65);
    legendMass->SetHeader("Legend");
    legendMass->AddEntry(PMass,"Parent Mass U+U","l");
    legendMass->AddEntry(PMassAu,"Parent Mass Au+Au","l");
    legendMass->Draw("same");
    gPad->Print("plots/plot_PMassBoth.png");


    makeCanvas2();
    TH1F *PMassRatio = (TH1F*)PMassAu->Clone("PMassRatio");
    PMassRatio->SetLineColor(kBlack);
    PMassRatio->Divide(PMass);
    PMassRatio->GetXaxis()->SetTitle("Pair Invariant Mass (GeV/c^{2})");
    PMassRatio->GetYaxis()->SetTitle("Pair Invariant Mass Ratio, Au+Au/U+U");
    PMassRatio->GetXaxis()->SetRangeUser(0.3,2);
    PMassRatio->GetYaxis()->SetRangeUser(0,2);
    PMassRatio->SetTitle("Parent Mass Ratio");
    PMassRatio->Draw();
    gPad->Print("plots/plot_PMassRatio.png");   


    makeCanvas2();
    mPt2->SetLineColor(kBlack);
    mPt2->GetXaxis()->SetTitle("Pt^{2} (GeV/c)^{2} U");
    mPt2->GetYaxis()->SetTitle("Counts (Normalized)");
    gPad->SetLogy();
    mPt2->Scale(1/mPt2->GetEntries());
    mPt2->Fit("expo", "R+", "", 0.001, .022);
    gStyle->SetOptFit(1111);
    mPt2->Draw();
    gPad->Print("plots/plot_mPt2.png");  

    makeCanvas2();
    mPt2Au->SetLineColor(kRed);
    mPt2Au->GetXaxis()->SetTitle("Pt^{2} (GeV/c)^{2} Au");
    mPt2Au->GetYaxis()->SetTitle("Counts (Normalized)");
    mPt2Au->Scale(1/mPt2Au->GetEntries());
    gPad->SetLogy();
    mPt2Au->Fit("expo", "R+", "", 0.001, .022);
    gStyle->SetOptFit(1111);
    mPt2Au->Draw();
    gPad->Print("plots/plot_mPt2Au.png");  

    makeCanvas2();
    TH1F *Pt2Ratio = (TH1F*)mPt2Au->Clone("Pt2Ratio");
    Pt2Ratio->SetLineColor(kBlack);
    Pt2Ratio->Divide(mPt2);
    Pt2Ratio->GetXaxis()->SetTitle("P_{T}^{2} (GeV/c^{2})^{2}");
    Pt2Ratio->GetYaxis()->SetTitle("P_{T}^{2} Ratio, Au+Au/U+U");
    Pt2Ratio->GetYaxis()->SetRangeUser(0,2);
    Pt2Ratio->SetTitle("Parent P_{T}^{2} Ratio");
    Pt2Ratio->Fit("pol1", "", "", 0, 0.012);
    Pt2Ratio->Draw();
    //auto * legendPt2 = new TLegend(0.79,0.5,.99,0.65);
    //legendPt2->SetHeader("Legend");
    //legendPt2->AddEntry(mPt2,"P_{T}^{2} U+U","l");
    //legendPt2->AddEntry(mPt2Au,"P_{T}^{2} Au+Au","l");
    //legendPt2->Draw("same");
    gPad->Print("plots/plot_Pt2Ratio.png");   


    makeCanvas2();
    mRapidity->SetLineColor(kBlack);
    mRapidity->Scale(1/mRapidity->GetEntries());
    mRapidity->Draw();

    mRapidityAu->SetLineColor(kRed);
    mRapidityAu->Scale(1/mRapidityAu->GetEntries());
    mRapidityAu->GetXaxis()->SetTitle("Rapidity");
    mRapidityAu->GetYaxis()->SetTitle("Counts");
    mRapidityAu->Draw("same");
    auto * legendY = new TLegend(0.79,0.5,.99,0.65);
    legendY->AddEntry(mRapidity, "Rapidity U+U", "l");
    legendY->AddEntry(mRapidityAu,"Rapidity Au+Au", "l");
    legendY->Draw("same");
    gPad->Print("plots/plot_rapidityBoth.png");


    makeCanvas2();
    TH1F *YRatio = (TH1F*)mRapidityAu->Clone("YRatio");
    YRatio->SetLineColor(kBlack);
    YRatio->Divide(mRapidity);
    YRatio->GetXaxis()->SetTitle("Rapidity");
    YRatio->GetYaxis()->SetTitle("Rapidity Ratio, Au+Au/U+U");
    YRatio->GetYaxis()->SetRangeUser(0,2);
    YRatio->SetTitle("Rapidity Ratio");
    //Pt2Ratio->Fit("pol1", "", "", 0, 0.012);
    YRatio->Draw();
    gPad->Print("plots/plot_rapidityRatio.png");


    makeCanvas2();
    m2Ptcos4phimomentsAu->SetTitle("scaled cos(4#phi) moments vs. P_{T} Au; P_{T} (GeV) Au; 2<cos(4#phi)>");
    m2Ptcos4phimomentsAu->GetXaxis()->SetTitle("<P_{T}> (GeV/c) Au");
    m2Ptcos4phimomentsAu->GetYaxis()->SetTitle("cos(4#phi) moments scaled");
    m2Ptcos4phimomentsAu->SetLineColor(kRed);
    m2Ptcos4phimomentsAu->Draw();
    //auto *phiFit4Au = new TF1("4phifit", "pol0");
    //m2Ptcos4phimomentsAu->Fit("4phifit", "", "", 0., 0.1);
    //gStyle->SetOptFit(1111);
    //phiFit4Au->Draw("same");
    //gPad->Print( "plots/plot_m2Ptcos4phimomentsAu.png");

    //makeCanvas2();
    m2Ptcos4phimoments->SetTitle("scaled cos(4#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(4#phi)>");
    m2Ptcos4phimoments->GetXaxis()->SetTitle("<P_{T}> (GeV/c)");
    m2Ptcos4phimoments->GetYaxis()->SetTitle("cos(4#phi) moments scaled");
    m2Ptcos4phimoments->SetLineColor(kBlack);
    m2Ptcos4phimoments->Draw("same");
    auto *phiFit4 = new TF1("4phifit", "pol0");
    //m2Ptcos4phimoments->Fit("4phifit", "", "", 0., 0.1);
    //gStyle->SetOptFit(1111);
    //phiFit4->Draw("same");
    auto * legendcos4phi = new TLegend(0.79,0.5,.99,0.65);
    legendcos4phi->AddEntry(m2Ptcos4phimomentsAu, "cos(4#phi) moments Au+Au", "l");
    legendcos4phi->AddEntry(m2Ptcos4phimoments,"cos(4#phi) moments U+U", "l");
    legendcos4phi->Draw("same");
    gPad->Print( "plots/plot_m2Ptcos4phimomentsBothWide.png");



    makeCanvas2();
    m2Ptcos3phimomentsAu->SetTitle("scaled cos(3#phi) moments vs. P_{T} Au; P_{T} (GeV) Au; 2<cos(3#phi)>");
    m2Ptcos3phimomentsAu->GetXaxis()->SetTitle("<P_{T}> (GeV/c) Au");
    m2Ptcos3phimomentsAu->GetYaxis()->SetTitle("cos(3#phi) moments scaled");
    m2Ptcos3phimomentsAu->SetLineColor(kBlack);
    m2Ptcos3phimomentsAu->Draw();
    auto *phiFit3Au = new TF1("3phifit", "pol0");
    m2Ptcos3phimomentsAu->Fit("3phifit", "", "", 0., 0.1);
    gStyle->SetOptFit(1111);
    phiFit3Au->Draw("same");
    gPad->Print( "plots/plot_m2Ptcos3phimomentsAuWide.png");


    makeCanvas2();
    m2Ptcos2phimomentsAu->SetTitle("scaled cos(2#phi) moments vs. P_{T} Au; P_{T} (GeV) Au; 2<cos(2#phi)>");
    m2Ptcos2phimomentsAu->GetXaxis()->SetTitle("<P_{T}> (GeV/c)");
    m2Ptcos2phimomentsAu->GetYaxis()->SetTitle("cos(2#phi) moments scaled");
    m2Ptcos2phimomentsAu->SetLineColor(kRed);
    m2Ptcos2phimomentsAu->Draw();
    //auto *phiFit2Au = new TF1("2phifit", "pol0");
    //m2Ptcos2phimomentsAu->Fit("2phifit", "", "", 0., 0.1);
    //gStyle->SetOptFit(1111);
    //phiFit2Au->Draw("same");
    //gPad->Print( "plots/plot_m2Ptcos2phimomentsAu.png");
     //makeCanvas2();
    m2Ptcos2phimoments->SetTitle("scaled cos(2#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
    //m2Ptcos2phimoments->GetXaxis()->SetTitle("<P_{T}> (GeV/c)");
    //m2Ptcos2phimoments->GetYaxis()->SetTitle("cos(2#phi) moments scaled");
    m2Ptcos2phimoments->SetLineColor(kBlack);
    m2Ptcos2phimoments->Draw("same");
    //auto *phiFit2 = new TF1("2phifit", "pol0");
    //m2Ptcos2phimoments->Fit("2phifit", "", "", 0., 0.1);
    //gStyle->SetOptFit(1111);
    //phiFit2->Draw("same");
    auto * legendcos2phi = new TLegend(0.79,0.5,.99,0.65);
    legendcos2phi->AddEntry(m2Ptcos2phimomentsAu, "cos(2#phi) moments Au+Au", "l");
    legendcos2phi->AddEntry(m2Ptcos2phimoments,"cos(2#phi) moments U+U", "l");
    legendcos2phi->Draw("same");
    gPad->Print( "plots/plot_m2Ptcos2phimomentsBothWide.png");



    makeCanvas2();
    m2PtcosphimomentsAu->SetTitle("scaled cos(#phi) moments vs. P_{T} Au; P_{T} (GeV) Au; 2<cos(#phi)>");
    m2PtcosphimomentsAu->GetXaxis()->SetTitle("<P_{T}> (GeV/c) Au");
    m2PtcosphimomentsAu->GetYaxis()->SetTitle("cos(#phi) moments scaled");
    m2PtcosphimomentsAu->SetLineColor(kBlack);
    m2PtcosphimomentsAu->Draw();
    auto *phiFitAu = new TF1("phifit", "pol0");
    m2PtcosphimomentsAu->Fit("phifit", "", "", 0., 0.1);
    gStyle->SetOptFit(1111);
    phiFitAu->Draw("same");
    gPad->Print( "plots/plot_m2PtcosphimomentsAuWide.png");


    makeCanvas2();
    TH2D *cos2phiPtmomentsRatio = (TH2D*)m2Ptcos2phimomentsAu->Clone("cos2phiPtmomentsRatio");
    cos2phiPtmomentsRatio->SetLineColor(kBlack);
    cos2phiPtmomentsRatio->Divide(m2Ptcos2phimoments);
    cos2phiPtmomentsRatio->GetXaxis()->SetTitle("<P_{T}> GeV/c");
    cos2phiPtmomentsRatio->GetYaxis()->SetTitle("Moment Ratio, Au+Au/U+U");
    cos2phiPtmomentsRatio->GetYaxis()->SetRangeUser(-3,3);
    cos2phiPtmomentsRatio->SetTitle("cos(2#phi) Moments Ratio");
    cos2phiPtmomentsRatio->Draw();
    gPad->Print("plots/plot_cos2phiPtmomentsRatioWide.png");

    makeCanvas2();
    TH2D *cos4phiPtmomentsRatio = (TH2D*)m2Ptcos4phimomentsAu->Clone("cos4phiPtmomentsRatio");
    cos4phiPtmomentsRatio->SetLineColor(kBlack);
    cos4phiPtmomentsRatio->Divide(m2Ptcos4phimoments);
    cos4phiPtmomentsRatio->GetXaxis()->SetTitle("<P_{T}> (GeV/c)");
    cos4phiPtmomentsRatio->GetYaxis()->SetTitle("Moment Ratio, Au+Au/U+U");
    cos4phiPtmomentsRatio->GetYaxis()->SetRangeUser(-4,4);
    cos4phiPtmomentsRatio->SetTitle("cos(4#phi) Moments Ratio");
    cos4phiPtmomentsRatio->Draw();
    gPad->Print("plots/plot_cos4phiPtmomentsRatioWide.png");

    

    makeCanvas2();
    m2Ptcos3phimoments->SetTitle("scaled cos(3#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(3#phi)>");
    m2Ptcos3phimoments->GetXaxis()->SetTitle("<P_{T}> (GeV/c)");
    m2Ptcos3phimoments->GetYaxis()->SetTitle("cos(3#phi) moments scaled");
    m2Ptcos3phimoments->SetLineColor(kBlack);
    m2Ptcos3phimoments->Draw();
    auto *phiFit3 = new TF1("3phifit", "pol0");
    m2Ptcos3phimoments->Fit("3phifit", "", "", 0., 0.1);
    gStyle->SetOptFit(1111);
    phiFit3->Draw("same");
    gPad->Print( "plots/plot_m2Ptcos3phimomentsWide.png");


    makeCanvas2();
    m2Ptcosphimoments->SetTitle("scaled cos(#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(#phi)>");
    m2Ptcosphimoments->GetXaxis()->SetTitle("<P_{T}> (GeV/c)");
    m2Ptcosphimoments->GetYaxis()->SetTitle("cos(#phi) moments scaled");
    m2Ptcosphimoments->SetLineColor(kBlack);
    m2Ptcosphimoments->Draw();
    auto *phiFit = new TF1("phifit", "pol0");
    m2Ptcosphimoments->Fit("phifit", "", "", 0., 0.1);
    gStyle->SetOptFit(1111);
    phiFit->Draw("same");
    gPad->Print( "plots/plot_m2PtcosphimomentsWide.png");



    auto *m2Ptcos4phimomentsHigherZDC = cos4phivPtHigherZDC->ProfileY("m2Ptcos4phimomentsHigherZDC", 1, -1);
    auto *m2Ptcos4phimomentsLowerZDC = cos4phivPtLowerZDC->ProfileY("m2Ptcos4phimomentsLowerZDC", 1, -1);
    makeCanvas2();
    m2Ptcos4phimomentsHigherZDC->SetLineColor(kBlack);
    m2Ptcos4phimomentsHigherZDC->Draw();
    m2Ptcos4phimomentsHigherZDC->GetXaxis()->SetTitle("<P_{T}> (GeV/c)");
    m2Ptcos4phimomentsHigherZDC->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsLowerZDC->SetLineColor(kRed);
    m2Ptcos4phimomentsLowerZDC->Draw("same");
    auto * legendcos4phiZDC = new TLegend(0.79,0.3,.99,0.55);
    legendcos4phiZDC->AddEntry(m2Ptcos4phimomentsHigherZDC, "cos(4#phi) moments \n U+U ZDC's > 200", "l");
    legendcos4phiZDC->AddEntry(m2Ptcos4phimomentsLowerZDC,"cos(4#phi) moments \n U+U ZDC's <= 200", "l");
    legendcos4phiZDC->SetTextSize(0.014);
    legendcos4phiZDC->Draw("same");
    gPad->Print("plots/plot_m2Ptcos4phimomentsZDCcutWide.png");

    
    auto *m2Ptcos4phimomentsHigherZDCAu = cos4phivPtHigherZDCAu->ProfileY("m2Ptcos4phimomentsHigherZDCAu", 1, -1);
    auto *m2Ptcos4phimomentsLowerZDCAu = cos4phivPtLowerZDCAu->ProfileY("m2Ptcos4phimomentsLowerZDCAu", 1, -1);
    makeCanvas2();
    m2Ptcos4phimomentsHigherZDCAu->SetLineColor(kBlack);
    m2Ptcos4phimomentsHigherZDCAu->Draw();
    m2Ptcos4phimomentsHigherZDCAu->GetXaxis()->SetTitle("<P_{T}> (GeV/c) Au");
    m2Ptcos4phimomentsHigherZDCAu->GetYaxis()->SetTitle("2<cos(4#phi)>");
    m2Ptcos4phimomentsLowerZDCAu->SetLineColor(kRed);
    m2Ptcos4phimomentsLowerZDCAu->Draw("same");
    auto * legendcos4phiZDCAu = new TLegend(0.79,0.3,.99,0.55);
    legendcos4phiZDCAu->AddEntry(m2Ptcos4phimomentsHigherZDCAu, "cos(4#phi) moments \n Au+Au ZDC's > 200", "l");
    legendcos4phiZDCAu->AddEntry(m2Ptcos4phimomentsLowerZDCAu,"cos(4#phi) moments \n Au+Au ZDC's <= 200", "l");
    legendcos4phiZDCAu->SetTextSize(0.014);
    legendcos4phiZDCAu->Draw("same");
    gPad->Print("plots/plot_m2Ptcos4phimomentsZDCcutAuWide.png");








}
