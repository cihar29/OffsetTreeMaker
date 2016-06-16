//Chad Harrington 3/23/2015

#include "TFile.h"
#include "TTree.h"
#include "plugins/parsePileUpJSON2.h"
using namespace std;

void fixmu(){

  parsePileUpJSON2("pileup_5_24_16.txt");
  
  TFile* file = TFile::Open("Data_tree.root");
  TTree* t = (TTree*) file->Get("T");

  t->SetBranchStatus("*",0);

  t->SetBranchStatus("run",1);
  t->SetBranchStatus("lumi",1);
  t->SetBranchStatus("bx",1);
  t->SetBranchStatus("event",1);
  t->SetBranchStatus("rho",1);
  t->SetBranchStatus("rhoC0",1);
  t->SetBranchStatus("rhoCC",1);
  t->SetBranchStatus("nPVall",1);
  t->SetBranchStatus("nPV",1);

  t->SetBranchStatus("pv_ndof",1);
  t->SetBranchStatus("pv_z",1);
  t->SetBranchStatus("pv_rho",1);
  t->SetBranchStatus("nEta",1);
  t->SetBranchStatus("energy",1);
  t->SetBranchStatus("et",1);
  t->SetBranchStatus("eRMS",1);

  t->SetBranchStatus("fchm",1);
  t->SetBranchStatus("fchu",1);
  t->SetBranchStatus("fnh",1);
  t->SetBranchStatus("fne",1);
  t->SetBranchStatus("fhfh",1);
  t->SetBranchStatus("fhfe",1);
  t->SetBranchStatus("flep",1);
  t->SetBranchStatus("funtrk",1);

  t->SetBranchStatus("ht",1);
  t->SetBranchStatus("nJets",1);
  t->SetBranchStatus("jet_eta",1);
  t->SetBranchStatus("jet_phi",1);
  t->SetBranchStatus("jet_pt",1);
  t->SetBranchStatus("jet_area",1);

  t->SetBranchStatus("jet_ch",1);
  t->SetBranchStatus("jet_nh",1);
  t->SetBranchStatus("jet_ne",1);
  t->SetBranchStatus("jet_hfh",1);
  t->SetBranchStatus("jet_hfe",1);
  t->SetBranchStatus("jet_lep",1);

  TFile *newfile = new TFile("muTree.root","recreate");
  TTree *newtree = t->CloneTree();
  Long64_t nEntries = newtree->GetEntries();

  float mu;
  TBranch* bmu = newtree->Branch("mu", &mu, "mu/F");
  int run, lumi;
  newtree->SetBranchAddress("run", &run);
  newtree->SetBranchAddress("lumi", &lumi);


  /*ULong64_t event;
  newtree->SetBranchAddress("event", &event);

  TFile* maskFile = TFile::Open("Data71300_tree.root");
  TTree* maskT = (TTree*) maskFile->Get("T");
  maskT->BuildIndex("run", "event");

  float tmu;
  maskT->SetBranchAddress("mu", &tmu);*/


  for (Long64_t n=0; n<nEntries; n++){
    newtree->GetEntry(n);

    mu = getAvgPU( run, lumi );

    //if (maskT->GetEntryWithIndex(run, event) == -1) {cout << run << "\t" << event << endl; continue;}
    //mu = tmu;

    bmu->Fill();
  }

  newtree->Print();
  newfile->Write();
  delete file;
  delete newfile;
}
