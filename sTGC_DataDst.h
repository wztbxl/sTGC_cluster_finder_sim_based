//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 26 20:07:54 2021 by ROOT version 6.18/00
// from TTree sTGC_DataDst/sTGC_DataDst
// found on file: test_Dignal_time_debug.root
//////////////////////////////////////////////////////////

#ifndef sTGC_DataDst_h
#define sTGC_DataDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class sTGC_DataDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           mEvtID;
   Int_t           nFireChannel;
   vector<int>     *Channel;
   vector<int>     *BCID;
   vector<int>     *ADC;

   // List of branches
   TBranch        *b_mEvtID;   //!
   TBranch        *b_nFireChannel;   //!
   TBranch        *b_Channel;   //!
   TBranch        *b_BCID;   //!
   TBranch        *b_ADC;   //!

   sTGC_DataDst(TTree *tree=0);
   virtual ~sTGC_DataDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef sTGC_DataDst_cxx
sTGC_DataDst::sTGC_DataDst(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../stgc-cluster-sim/test_Dignal_time_debug.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test_Dignal_time_debug.root");
      }
      f->GetObject("sTGC_DataDst",tree);

   }
   Init(tree);
}

sTGC_DataDst::~sTGC_DataDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sTGC_DataDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sTGC_DataDst::LoadTree(Long64_t entry)
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

void sTGC_DataDst::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Channel = 0;
   BCID = 0;
   ADC = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mEvtID", &mEvtID, &b_mEvtID);
   fChain->SetBranchAddress("nFireChannel", &nFireChannel, &b_nFireChannel);
   fChain->SetBranchAddress("Channel", &Channel, &b_Channel);
   fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
   fChain->SetBranchAddress("ADC", &ADC, &b_ADC);
   Notify();
}

Bool_t sTGC_DataDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sTGC_DataDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sTGC_DataDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sTGC_DataDst_cxx
