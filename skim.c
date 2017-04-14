//Chad Harrington 3/23/2015
//Skim a TTree

#include "TFile.h"
#include "TTree.h"
using namespace std;

void skim(){
  
  int numSkip = 1013;

  TFile* file = TFile::Open("bigFile.root");

  TTree* tree = (TTree*) file->Get("T");
  Long64_t nEntries = tree->GetEntries();

  ULong64_t event;
  tree->SetBranchAddress("event", &event);
  tree->SetBranchStatus("*",1);

  TFile *newfile = new TFile("skimmedTree.root","recreate");
  TTree *newtree = tree->CloneTree(0);

  for (Long64_t n=0; n<nEntries; n++) {
    tree->GetEntry(n);

    if (event % numSkip == 0) newtree->Fill();
  }

  newtree->Print();
  newfile->Write();
  delete file;
  delete newfile;
}
