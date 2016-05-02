//Chad Harrington 3/23/2015
//MC Pileup reweighting

#include "TFile.h"
#include "TH1.h"
using namespace std;

void weights(){
  
  TFile* mcFile = TFile::Open("MC76x_mu.root");
  TFile* dataFile = TFile::Open("Data76x_mu.root");

  TString hname;

  hname = "nPU";

  TH1F* h_MC = (TH1F*) mcFile->Get(hname);
  TH1F* h_Data = (TH1F*) dataFile->Get(hname);

  TH1F* h_MC_Data = (TH1F*) h_Data->Clone("h_MC_Data");
  h_MC_Data->Divide(h_MC);
  h_MC_Data->Scale( 1/ h_MC_Data->GetMaximum() );

  h_Data->Scale( 1 / h_Data->GetMaximum() );
  h_MC->Scale( 1 / h_MC->GetMaximum() );

  h_Data->Draw();
  h_MC->Draw("same");
  h_MC_Data->Draw("same");

  TString array = "      {";
  int size = h_MC_Data->GetNbinsX();

  for (int i=1; i<=size; i++){
  
    if ( i != 1 && i != size && (i-1)%10 == 0)
      array += "\n       ";

    array += Form( "%1.5f, ", h_MC_Data->GetBinContent(i) );
  }

  array += "}";
  cout << array << endl;
}
