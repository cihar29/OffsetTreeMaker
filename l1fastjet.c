//Chad Harrington 7/28/2015, root -l -b -q l1fastjet.c

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

void setStyle();

void l1fastjet(TString pf_type="all"){

  const int ETA_BINS = 82;
  double etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

  int Rlabel = 4;
  bool nPU_derived = true;
  bool rhoCentral = false;

  TFile* data_root = TFile::Open( Form("Data7716_R%i.root", Rlabel) );
  TFile* mc_root = TFile::Open( Form("MC7716_R%i.root", Rlabel) );

  ifstream scale_file( Form("plots/scalefactor/Spring16_25nsV1_DataMcSF_L1RC_AK%iPF", Rlabel) + pf_type + ".txt" );

  ifstream mc_file( Form("Spring16_25nsV1_MC_L1FastJet_AK%iPF", Rlabel) + pf_type + ".txt" );
  string scale_line, mc_line;

  //read first line
  getline(scale_file, scale_line);
  getline(mc_file, mc_line);

  double scale_p0[ETA_BINS], scale_p1[ETA_BINS], scale_p2[ETA_BINS];
  double mc_p0[ETA_BINS], mc_p1[ETA_BINS];

  for (int i=0; getline(scale_file,scale_line); i++){

    string str;
    int delim_pos;

    while (scale_line.at(0) == ' ')  //eta values have space in front of them
      scale_line.erase(0, 1);

    //loop over strings in data line
    for (int col_num=0; (delim_pos = scale_line.find(' ')) != -1; col_num++){

      str = scale_line.substr(0, delim_pos);
      scale_line.erase(0, delim_pos + 1);

      while (scale_line.at(0) == ' ')  //get rid of white space between columns
        scale_line.erase(0, 1);

      if (col_num == 5) scale_p0[i] = stod(str);
      else if (col_num == 6) scale_p1[i] = stod(str);
    }
    //last column after loop
    scale_p2[i] = stod(scale_line);
  }
  scale_file.close();

  int lines;
  for (int i=0; getline(mc_file,mc_line); i++){

    lines = i;
    string str;
    int delim_pos;

    while (mc_line.at(0) == ' ')  //eta values have space in front of them
      mc_line.erase(0, 1);

    //loop over strings in mc line
    for (int col_num=0; (delim_pos = mc_line.find(' ')) != -1; col_num++){

      str = mc_line.substr(0, delim_pos);
      mc_line.erase(0, delim_pos + 1);

      while (!mc_line.empty() && mc_line.at(0) == ' ')  //get rid of white space between columns
        mc_line.erase(0, 1);

      if (col_num == 9) mc_p0[i] = stod(str);
      else if (col_num == 10) mc_p1[i] = stod(str);
    }
  }

  TString hname;

  if (nPU_derived) hname = "nPU";
  else hname = "nPV";

  TH1F* h_nPU = (TH1F*) data_root->Get(hname);
  double mean_nPU = h_nPU->GetMean();

  int mean_bin = h_nPU->FindBin(mean_nPU);
  int low_bin, high_bin;
  double total_area = h_nPU->Integral();

  for (low_bin = mean_bin; (h_nPU->Integral(low_bin, mean_bin) / total_area) < 0.34; low_bin--){}
  for (high_bin = mean_bin; (h_nPU->Integral(mean_bin, high_bin) / total_area) < 0.34; high_bin++){}

  double nPU_low = h_nPU->GetBinCenter(low_bin);
  double nPU_high = h_nPU->GetBinCenter(high_bin);

  if (nPU_derived){
      if (rhoCentral) hname = "p_rhoCentral0_nPU";
      else hname = "p_rho_nPU";
  }
  else{
      if (rhoCentral) hname = "p_rhoCentral0_nPV";
      else hname = "p_rho_nPV";
  }

  TProfile* data_rho_nPU = (TProfile*) data_root->Get(hname);
  TProfile* mc_rho_nPU = (TProfile*) mc_root->Get(hname);

  double rho_nominal = data_rho_nPU->GetBinContent( data_rho_nPU->FindBin( mean_nPU ) );
  double rho_low = data_rho_nPU->GetBinContent( data_rho_nPU->FindBin( nPU_low ) );
  double rho_high = data_rho_nPU->GetBinContent( data_rho_nPU->FindBin( nPU_high ) );

  double mc_rho_mean = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin( mean_nPU ) );
  double mc_rho_low = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin( nPU_low ) );
  double mc_rho_high = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin( nPU_high ) );

  double eta[ETA_BINS], sf[ETA_BINS], sf_low[ETA_BINS], sf_high[ETA_BINS], new_p0[ETA_BINS], new_p1[ETA_BINS];

  for (int i=0; i<ETA_BINS; i++){

    eta[i] = 0.5*(etabins[i]+etabins[i+1]);

    sf_low[i] = scale_p0[i] + scale_p1[i]*rho_low + scale_p2[i]*rho_low*rho_low;
    sf_high[i] = scale_p0[i] + scale_p1[i]*rho_high + scale_p2[i]*rho_high*rho_high;

    sf[i] = scale_p0[i] + scale_p1[i]*rho_nominal + scale_p2[i]*rho_nominal*rho_nominal;
    new_p0[i] = sf[i] * mc_p0[i];
    new_p1[i] = sf[i] * mc_p1[i] * mc_rho_mean / rho_nominal;
  }

  ofstream writeFile( Form("Spring16_25nsV1_DATA_L1FastJet_AK%iPF", Rlabel) + pf_type + ".txt" );
  mc_file.clear();
  mc_file.seekg(0, mc_file.beg);

  //write first line
  getline(mc_file,mc_line);
  writeFile << mc_line << endl;

  for (int i=0; getline(mc_file,mc_line); i++){

    string str;
    int delim_pos;

    while (mc_line.at(0) == ' ')  //eta values have space in front of them
      mc_line.erase(0, 1);

    //loop over strings in mc line
    for (int col_num=0; (delim_pos = mc_line.find(' ')) != -1; col_num++){

      str = mc_line.substr(0, delim_pos);
      mc_line.erase(0, delim_pos + 1);

      while (!mc_line.empty() && mc_line.at(0) == ' ')  //get rid of white space between columns
        mc_line.erase(0, 1);

      if (col_num == 9) str = to_string(new_p0[i]);
      else if (col_num == 10) str = to_string(new_p1[i]);

      writeFile << str << setw(15);
    }
    writeFile << mc_line << endl;
  }
  mc_file.close();
  writeFile.close();

  setStyle();

  TGraphErrors* graph = new TGraphErrors(ETA_BINS, eta, sf, 0, 0);
  TGraphErrors* g_high = new TGraphErrors(ETA_BINS, eta, sf_high, 0, 0);
  TGraphErrors* g_low = new TGraphErrors(ETA_BINS, eta, sf_low, 0, 0);  

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  TH1F* h = new TH1F("h", "h", ETA_BINS, etabins);

  h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("Scale Factor");
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetRangeUser(0.9, 1.3);
  h->Draw();

  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(1);
  graph->Draw("psame");

  g_high->SetMarkerStyle(2);
  g_high->SetMarkerColor(kRed);
  g_high->Draw("psame");

  g_low->SetMarkerStyle(5);
  g_low->SetMarkerColor(kBlue);
  g_low->Draw("psame");

  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.04);

  if (pf_type.EqualTo("all"))
    text.DrawLatex(0.17, 0.96, Form("AK%i PF", Rlabel));
  else
    text.DrawLatex(0.17, 0.96, Form("AK%i PF%s", Rlabel, pf_type.Data()));

  string var_string;
  if (nPU_derived) var_string = "#mu";
  else var_string = "N_{PV}";

  string rho_string = "";
  if (rhoCentral) rho_string = "central, ";

  float leg = 0.35;

  text.DrawLatex(0.45, leg, Form( "#bf{%s = %4.1f^{ #color[2]{%4.1f}}_{ #color[4]{%4.1f}}}" , var_string.c_str(), mean_nPU, nPU_high, nPU_low ) );
  text.DrawLatex(0.45, leg-0.07, Form( "#bf{#rho_{%sData} = %4.1f^{ #color[2]{%4.1f}}_{ #color[4]{%4.1f}} GeV}" , rho_string.c_str(), rho_nominal, rho_high, rho_low ) );
  text.DrawLatex(0.45, leg-0.14, Form( "#bf{#rho_{%sMC} = %4.1f^{ #color[2]{%4.1f}}_{ #color[4]{%4.1f}} GeV}" , rho_string.c_str(), mc_rho_mean, mc_rho_high, mc_rho_low ) );

  text.SetTextSize(0.035);
  text.SetTextColor(1);
  text.SetTextFont(42);
  text.DrawLatex(0.53, 0.96, "Run 2016B - 6.26 fb^{-1} (13 TeV)");

  c->Print("l1fastjet_scalefactor_eta_PF" + pf_type + ".pdf");
}

void setStyle(){

//Style//

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.9);

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  //tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->SetOptStat(0);

  tdrStyle->cd();

//End Style//
}
