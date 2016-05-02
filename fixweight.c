#include "TFile.h"
#include "TTree.h"

void fixweight(){

  TFile* oldfile = TFile::Open("temp.root");  
  TTree* oldtree = (TTree*) oldfile->Get("T");

  oldtree->SetBranchStatus("*",0);
  oldtree->SetBranchStatus("mu",1);
  oldtree->SetBranchStatus("rho",1);
  oldtree->SetBranchStatus("rhoC0",1);
  oldtree->SetBranchStatus("rhoCC",1);
  oldtree->SetBranchStatus("nPVall",1);
  oldtree->SetBranchStatus("nPV",1);

  oldtree->SetBranchStatus("pv_ndof",1);
  oldtree->SetBranchStatus("pv_z",1);
  oldtree->SetBranchStatus("pv_rho",1);
  oldtree->SetBranchStatus("nEta",1);
  oldtree->SetBranchStatus("energy",1);
  oldtree->SetBranchStatus("et",1);
  oldtree->SetBranchStatus("eRMS",1);

  oldtree->SetBranchStatus("fchm",1);
  oldtree->SetBranchStatus("fchu",1);
  oldtree->SetBranchStatus("fnh",1);
  oldtree->SetBranchStatus("fne",1);
  oldtree->SetBranchStatus("fhfh",1);
  oldtree->SetBranchStatus("fhfe",1);
  oldtree->SetBranchStatus("flep",1);
  oldtree->SetBranchStatus("funtrk",1);

  oldtree->SetBranchStatus("ht",1);
  oldtree->SetBranchStatus("nJets",1);
  oldtree->SetBranchStatus("jet_eta",1);
  oldtree->SetBranchStatus("jet_phi",1);
  oldtree->SetBranchStatus("jet_pt",1);
  oldtree->SetBranchStatus("jet_area",1);

  oldtree->SetBranchStatus("jet_ch",1);
  oldtree->SetBranchStatus("jet_nh",1);
  oldtree->SetBranchStatus("jet_ne",1);
  oldtree->SetBranchStatus("jet_hfh",1);
  oldtree->SetBranchStatus("jet_hfe",1);
  oldtree->SetBranchStatus("jet_lep",1);

  TFile *newfile = new TFile("MC101_76x.root","recreate");
  TTree *newtree = oldtree->CloneTree();
  Long64_t nEntries = newtree->GetEntries();

  float mu, weight;
  TBranch* bwgt = newtree->Branch("weight", &weight, "weight/F");
  newtree->SetBranchAddress("mu", &mu);

  TH1F* h = new TH1F("mu", "mu", 100, 0, 50);
  float weights[] =
      {1.00000, 0.22358, 0.08857, 0.10841, 0.11072, 0.05502, 0.07544, 0.03362, 0.07415, 0.03482, 
       0.01966, 0.02403, 0.00540, 0.00133, 0.00053, 0.00089, 0.01576, 0.02411, 0.02743, 0.03908, 
       0.02668, 0.03015, 0.03514, 0.04335, 0.03932, 0.02878, 0.03152, 0.01516, 0.00320, 0.00043, 
       0.00000, 0.00002, 0.00003, 0.00003, 0.00260, 0.00419, 0.00873, 0.00742, 0.00670, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 
       0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000 };

  for (Long64_t n=0; n<nEntries; n++){
    newtree->GetEntry(n);

    int index = h->FindBin(mu);
    weight = weights[index-1];

    bwgt->Fill();
  }
  delete h;

  newtree->Print();
  newfile->Write();
  delete oldfile;
  delete newfile;
}
