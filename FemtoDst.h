//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 15 10:36:07 2023 by ROOT version 6.28/00
// from TTree FemtoDst/FemtoDst with MC info
// found on file: /Users/Nick/Downloads/FemtoDst_Run12UU_wZDC.root
//////////////////////////////////////////////////////////

#ifndef FemtoDst_h
#define FemtoDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"
#include "TClonesArray.h"

class FemtoDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxTracks = 2;
   static constexpr Int_t kMaxBTofPidTraits = 2;

   // Declaration of leaf types
 //FemtoEvent      *Event;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Float_t         mPrimaryVertex_mX1;
   Float_t         mPrimaryVertex_mX2;
   Float_t         mPrimaryVertex_mX3;
   Float_t         mVpdVz;
   Int_t           mRunId;
   Int_t           mEventId;
   UShort_t        mGRefMult;
   UShort_t        mRefMult;
   UShort_t        mNBTofMatched;
   UShort_t        mNBTofHits;
   UShort_t        mNTracks;
   UShort_t        mNVerts;
   Float_t         mRanking;
   UShort_t        mNVpdEast;
   UShort_t        mNVpdWest;
   UShort_t        mRunIndex;
   UShort_t        mZDCEast;
   UShort_t        mZDCWest;
   Int_t           Tracks_;
   UInt_t          Tracks_fUniqueID[kMaxTracks];   //[Tracks_]
   UInt_t          Tracks_fBits[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mPt[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mEta[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mPhi[kMaxTracks];   //[Tracks_]
   UShort_t        Tracks_mId[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mDedx[kMaxTracks];   //[Tracks_]
   Char_t          Tracks_mQ[kMaxTracks];   //[Tracks_]
   Char_t          Tracks_mNHitsFit[kMaxTracks];   //[Tracks_]
   Char_t          Tracks_mNHitsMax[kMaxTracks];   //[Tracks_]
   UChar_t         Tracks_mNHitsDedx[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mNSigmaPion[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mNSigmaKaon[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mNSigmaProton[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mNSigmaElectron[kMaxTracks];   //[Tracks_]
   Float_t         Tracks_mDCA[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mBTofPidTraitsIndex[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mMtdPidTraitsIndex[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mEmcPidTraitsIndex[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mHelixIndex[kMaxTracks];   //[Tracks_]
   Short_t         Tracks_mMcIndex[kMaxTracks];   //[Tracks_]
   Int_t           BTofPidTraits_;
   UInt_t          BTofPidTraits_fUniqueID[kMaxBTofPidTraits];   //[BTofPidTraits_]
   UInt_t          BTofPidTraits_fBits[kMaxBTofPidTraits];   //[BTofPidTraits_]
   Float_t         BTofPidTraits_mBTofBeta[kMaxBTofPidTraits];   //[BTofPidTraits_]
   Float_t         BTofPidTraits_mBTof[kMaxBTofPidTraits];   //[BTofPidTraits_]
   Float_t         BTofPidTraits_mLength[kMaxBTofPidTraits];   //[BTofPidTraits_]
   Float_t         BTofPidTraits_mBTofYLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
   Float_t         BTofPidTraits_mBTofZLocal[kMaxBTofPidTraits];   //[BTofPidTraits_]
   UChar_t         BTofPidTraits_mBTofMatchFlag[kMaxBTofPidTraits];   //[BTofPidTraits_]
   Short_t         BTofPidTraits_mIdTruth[kMaxBTofPidTraits];   //[BTofPidTraits_]

   // List of branches
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_mPrimaryVertex_mX1;   //!
   TBranch        *b_Event_mPrimaryVertex_mX2;   //!
   TBranch        *b_Event_mPrimaryVertex_mX3;   //!
   TBranch        *b_Event_mVpdVz;   //!
   TBranch        *b_Event_mRunId;   //!
   TBranch        *b_Event_mEventId;   //!
   TBranch        *b_Event_mGRefMult;   //!
   TBranch        *b_Event_mRefMult;   //!
   TBranch        *b_Event_mNBTofMatched;   //!
   TBranch        *b_Event_mNBTofHits;   //!
   TBranch        *b_Event_mNTracks;   //!
   TBranch        *b_Event_mNVerts;   //!
   TBranch        *b_Event_mRanking;   //!
   TBranch        *b_Event_mNVpdEast;   //!
   TBranch        *b_Event_mNVpdWest;   //!
   TBranch        *b_Event_mRunIndex;   //!
   TBranch        *b_Event_mZDCEast;   //!
   TBranch        *b_Event_mZDCWest;   //!
   TBranch        *b_Tracks_;   //!
   TBranch        *b_Tracks_fUniqueID;   //!
   TBranch        *b_Tracks_fBits;   //!
   TBranch        *b_Tracks_mPt;   //!
   TBranch        *b_Tracks_mEta;   //!
   TBranch        *b_Tracks_mPhi;   //!
   TBranch        *b_Tracks_mId;   //!
   TBranch        *b_Tracks_mDedx;   //!
   TBranch        *b_Tracks_mQ;   //!
   TBranch        *b_Tracks_mNHitsFit;   //!
   TBranch        *b_Tracks_mNHitsMax;   //!
   TBranch        *b_Tracks_mNHitsDedx;   //!
   TBranch        *b_Tracks_mNSigmaPion;   //!
   TBranch        *b_Tracks_mNSigmaKaon;   //!
   TBranch        *b_Tracks_mNSigmaProton;   //!
   TBranch        *b_Tracks_mNSigmaElectron;   //!
   TBranch        *b_Tracks_mDCA;   //!
   TBranch        *b_Tracks_mBTofPidTraitsIndex;   //!
   TBranch        *b_Tracks_mMtdPidTraitsIndex;   //!
   TBranch        *b_Tracks_mEmcPidTraitsIndex;   //!
   TBranch        *b_Tracks_mHelixIndex;   //!
   TBranch        *b_Tracks_mMcIndex;   //!
   TBranch        *b_BTofPidTraits_;   //!
   TBranch        *b_BTofPidTraits_fUniqueID;   //!
   TBranch        *b_BTofPidTraits_fBits;   //!
   TBranch        *b_BTofPidTraits_mBTofBeta;   //!
   TBranch        *b_BTofPidTraits_mBTof;   //!
   TBranch        *b_BTofPidTraits_mLength;   //!
   TBranch        *b_BTofPidTraits_mBTofYLocal;   //!
   TBranch        *b_BTofPidTraits_mBTofZLocal;   //!
   TBranch        *b_BTofPidTraits_mBTofMatchFlag;   //!
   TBranch        *b_BTofPidTraits_mIdTruth;   //!

   FemtoDst(TTree *tree=0);
   virtual ~FemtoDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FemtoDst_cxx
FemtoDst::FemtoDst(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/Nick/Downloads/FemtoDst_Run12UU_wZDC.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Users/Nick/Downloads/FemtoDst_Run12UU_wZDC.root");
      }
      f->GetObject("FemtoDst",tree);

   }
   Init(tree);
}

FemtoDst::~FemtoDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FemtoDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FemtoDst::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FemtoDst::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Event_fBits);
   fChain->SetBranchAddress("mPrimaryVertex_mX1", &mPrimaryVertex_mX1, &b_Event_mPrimaryVertex_mX1);
   fChain->SetBranchAddress("mPrimaryVertex_mX2", &mPrimaryVertex_mX2, &b_Event_mPrimaryVertex_mX2);
   fChain->SetBranchAddress("mPrimaryVertex_mX3", &mPrimaryVertex_mX3, &b_Event_mPrimaryVertex_mX3);
   fChain->SetBranchAddress("mVpdVz", &mVpdVz, &b_Event_mVpdVz);
   fChain->SetBranchAddress("mRunId", &mRunId, &b_Event_mRunId);
   fChain->SetBranchAddress("mEventId", &mEventId, &b_Event_mEventId);
   fChain->SetBranchAddress("mGRefMult", &mGRefMult, &b_Event_mGRefMult);
   fChain->SetBranchAddress("mRefMult", &mRefMult, &b_Event_mRefMult);
   fChain->SetBranchAddress("mNBTofMatched", &mNBTofMatched, &b_Event_mNBTofMatched);
   fChain->SetBranchAddress("mNBTofHits", &mNBTofHits, &b_Event_mNBTofHits);
   fChain->SetBranchAddress("mNTracks", &mNTracks, &b_Event_mNTracks);
   fChain->SetBranchAddress("mNVerts", &mNVerts, &b_Event_mNVerts);
   fChain->SetBranchAddress("mRanking", &mRanking, &b_Event_mRanking);
   fChain->SetBranchAddress("mNVpdEast", &mNVpdEast, &b_Event_mNVpdEast);
   fChain->SetBranchAddress("mNVpdWest", &mNVpdWest, &b_Event_mNVpdWest);
   fChain->SetBranchAddress("mRunIndex", &mRunIndex, &b_Event_mRunIndex);
   fChain->SetBranchAddress("mZDCEast", &mZDCEast, &b_Event_mZDCEast);
   fChain->SetBranchAddress("mZDCWest", &mZDCWest, &b_Event_mZDCWest);
   fChain->SetBranchAddress("Tracks", &Tracks_, &b_Tracks_);
   fChain->SetBranchAddress("Tracks.fUniqueID", Tracks_fUniqueID, &b_Tracks_fUniqueID);
   fChain->SetBranchAddress("Tracks.fBits", Tracks_fBits, &b_Tracks_fBits);
   fChain->SetBranchAddress("Tracks.mPt", Tracks_mPt, &b_Tracks_mPt);
   fChain->SetBranchAddress("Tracks.mEta", Tracks_mEta, &b_Tracks_mEta);
   fChain->SetBranchAddress("Tracks.mPhi", Tracks_mPhi, &b_Tracks_mPhi);
   fChain->SetBranchAddress("Tracks.mId", Tracks_mId, &b_Tracks_mId);
   fChain->SetBranchAddress("Tracks.mDedx", Tracks_mDedx, &b_Tracks_mDedx);
   fChain->SetBranchAddress("Tracks.mQ", Tracks_mQ, &b_Tracks_mQ);
   fChain->SetBranchAddress("Tracks.mNHitsFit", Tracks_mNHitsFit, &b_Tracks_mNHitsFit);
   fChain->SetBranchAddress("Tracks.mNHitsMax", Tracks_mNHitsMax, &b_Tracks_mNHitsMax);
   fChain->SetBranchAddress("Tracks.mNHitsDedx", Tracks_mNHitsDedx, &b_Tracks_mNHitsDedx);
   fChain->SetBranchAddress("Tracks.mNSigmaPion", Tracks_mNSigmaPion, &b_Tracks_mNSigmaPion);
   fChain->SetBranchAddress("Tracks.mNSigmaKaon", Tracks_mNSigmaKaon, &b_Tracks_mNSigmaKaon);
   fChain->SetBranchAddress("Tracks.mNSigmaProton", Tracks_mNSigmaProton, &b_Tracks_mNSigmaProton);
   fChain->SetBranchAddress("Tracks.mNSigmaElectron", Tracks_mNSigmaElectron, &b_Tracks_mNSigmaElectron);
   fChain->SetBranchAddress("Tracks.mDCA", Tracks_mDCA, &b_Tracks_mDCA);
   fChain->SetBranchAddress("Tracks.mBTofPidTraitsIndex", Tracks_mBTofPidTraitsIndex, &b_Tracks_mBTofPidTraitsIndex);
   fChain->SetBranchAddress("Tracks.mMtdPidTraitsIndex", Tracks_mMtdPidTraitsIndex, &b_Tracks_mMtdPidTraitsIndex);
   fChain->SetBranchAddress("Tracks.mEmcPidTraitsIndex", Tracks_mEmcPidTraitsIndex, &b_Tracks_mEmcPidTraitsIndex);
   fChain->SetBranchAddress("Tracks.mHelixIndex", Tracks_mHelixIndex, &b_Tracks_mHelixIndex);
   fChain->SetBranchAddress("Tracks.mMcIndex", Tracks_mMcIndex, &b_Tracks_mMcIndex);
   fChain->SetBranchAddress("BTofPidTraits", &BTofPidTraits_, &b_BTofPidTraits_);
   fChain->SetBranchAddress("BTofPidTraits.fUniqueID", BTofPidTraits_fUniqueID, &b_BTofPidTraits_fUniqueID);
   fChain->SetBranchAddress("BTofPidTraits.fBits", BTofPidTraits_fBits, &b_BTofPidTraits_fBits);
   fChain->SetBranchAddress("BTofPidTraits.mBTofBeta", BTofPidTraits_mBTofBeta, &b_BTofPidTraits_mBTofBeta);
   fChain->SetBranchAddress("BTofPidTraits.mBTof", BTofPidTraits_mBTof, &b_BTofPidTraits_mBTof);
   fChain->SetBranchAddress("BTofPidTraits.mLength", BTofPidTraits_mLength, &b_BTofPidTraits_mLength);
   fChain->SetBranchAddress("BTofPidTraits.mBTofYLocal", BTofPidTraits_mBTofYLocal, &b_BTofPidTraits_mBTofYLocal);
   fChain->SetBranchAddress("BTofPidTraits.mBTofZLocal", BTofPidTraits_mBTofZLocal, &b_BTofPidTraits_mBTofZLocal);
   fChain->SetBranchAddress("BTofPidTraits.mBTofMatchFlag", BTofPidTraits_mBTofMatchFlag, &b_BTofPidTraits_mBTofMatchFlag);
   fChain->SetBranchAddress("BTofPidTraits.mIdTruth", BTofPidTraits_mIdTruth, &b_BTofPidTraits_mIdTruth);
   Notify();
}

Bool_t FemtoDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FemtoDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FemtoDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FemtoDst_cxx
