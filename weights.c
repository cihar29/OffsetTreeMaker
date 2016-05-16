//Chad Harrington 3/23/2015
//MC Pileup reweighting

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
using namespace std;

void weights(){
  
  TFile* mcFile = TFile::Open("MCTree.root");
  TFile* dataFile = TFile::Open("DataTree.root");

  TH1F* h_MC = new TH1F("h_mc","h_mc",100,0,50);
  TH1F* h_Data = new TH1F("h_data","h_data",100,0,50);
  float mcMu, dataMu;

  TTree* mcT = (TTree*) mcFile->Get("T");
  Long64_t mcEnt = mcT->GetEntries();
  cout << "Filling mu histograms...\nMC: " << mcEnt << " events" << endl;
  mcT->SetBranchAddress("mu", &mcMu);

  for (Long64_t n=0; n<mcEnt; n++){
    mcT->GetEntry(n);
    h_MC->Fill(mcMu);
  }

  TTree* dataT = (TTree*) dataFile->Get("T");
  Long64_t dataEnt = dataT->GetEntries();
  cout << "Data: " << dataEnt << " events" << endl;
  dataT->SetBranchAddress("mu", &dataMu);

  for (Long64_t n=0; n<dataEnt; n++){
    dataT->GetEntry(n);
    h_Data->Fill(dataMu);
  }

  cout << "Done\nFilling branches..." << endl;

  mcT->SetBranchStatus("*",0);
  mcT->SetBranchStatus("mu",1);
  mcT->SetBranchStatus("rho",1);
  mcT->SetBranchStatus("rhoC0",1);
  mcT->SetBranchStatus("rhoCC",1);
  mcT->SetBranchStatus("nPVall",1);
  mcT->SetBranchStatus("nPV",1);

  mcT->SetBranchStatus("pv_ndof",1);
  mcT->SetBranchStatus("pv_z",1);
  mcT->SetBranchStatus("pv_rho",1);
  mcT->SetBranchStatus("nEta",1);
  mcT->SetBranchStatus("energy",1);
  mcT->SetBranchStatus("et",1);
  mcT->SetBranchStatus("eRMS",1);

  mcT->SetBranchStatus("fchm",1);
  mcT->SetBranchStatus("fchu",1);
  mcT->SetBranchStatus("fnh",1);
  mcT->SetBranchStatus("fne",1);
  mcT->SetBranchStatus("fhfh",1);
  mcT->SetBranchStatus("fhfe",1);
  mcT->SetBranchStatus("flep",1);
  mcT->SetBranchStatus("funtrk",1);

  mcT->SetBranchStatus("ht",1);
  mcT->SetBranchStatus("nJets",1);
  mcT->SetBranchStatus("jet_eta",1);
  mcT->SetBranchStatus("jet_phi",1);
  mcT->SetBranchStatus("jet_pt",1);
  mcT->SetBranchStatus("jet_area",1);

  mcT->SetBranchStatus("jet_ch",1);
  mcT->SetBranchStatus("jet_nh",1);
  mcT->SetBranchStatus("jet_ne",1);
  mcT->SetBranchStatus("jet_hfh",1);
  mcT->SetBranchStatus("jet_hfe",1);
  mcT->SetBranchStatus("jet_lep",1);

  TFile *newfile = new TFile("newTree.root","recreate");
  TTree *newtree = mcT->CloneTree();
  Long64_t nEntries = newtree->GetEntries();

  float mu, weight;
  TBranch* bwgt = newtree->Branch("weight", &weight, "weight/F");
  newtree->SetBranchAddress("mu", &mu);

  h_Data->Divide(h_MC);
  h_Data->Scale( 1/ h_Data->GetMaximum() );

  for (Long64_t n=0; n<nEntries; n++){
    newtree->GetEntry(n);

    weight = float( h_Data->GetBinContent( h_Data->FindBin(mu) ) );
    bwgt->Fill();
  }
  delete dataFile, h_Data, h_MC;

  newtree->Print();
  newfile->Write();
  delete mcFile;
  delete newfile;
}
