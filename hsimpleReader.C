/// \file
/// \ingroup tutorial_tree
/// \notebook
/// TTreeReader simplest example.
///
/// Read data from hsimple.root (written by hsimple.C)
///
/// \macro_code
///
/// \author Anders Eie, 2013

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void hsimpleReader() {
   // Create a histogram for the values we read.
   auto myHist = new TH1F("h1","Vr",100,0,2);

   // Open the file containing the tree.
   auto myFile = TFile::Open("/Users/Nick/Desktop/Spring2023/FemtoDst_Run12UU_wZDC.root");
   if (!myFile || myFile->IsZombie()) {
      return;
   }
   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in.
   TTreeReader myReader("FemtoDst", myFile);

   // The branch "px" contains floats; access them as myPx.
   TTreeReaderValue<Float_t> myVr(myReader, "mVr");
   // The branch "py" contains floats, too; access those as myPy.
   //TTreeReaderValue<Float_t> myPy(myReader, "py");

   // Loop over all entries of the TTree or TChain.
   while (myReader.Next()) {
      // Just access the data as if myPx and myPy were iterators (note the '*'
      // in front of them):
      myHist->Fill(*myVr);
   }

   myHist->Draw();
}
