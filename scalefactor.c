//Chad Harrington 7/28/2015, root -l -b -q scalefactor.c

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

void setStyle();

void scalefactor(TString pf_type="all", int n1=1, int n2=25){

  setStyle();

  int Rlabel = 4;
  TFile* data_root = TFile::Open( Form("Data_R%i.root", Rlabel) );
  TFile* mc_root = TFile::Open( Form("MC_R%i.root", Rlabel) );

  bool nPU_derived = true;
  bool rhoCentral = false;

  ifstream data_file("./plots/indirectRho/" + pf_type + Form("/Spring16_25nsV1_DATA_L1RC_AK%iPF", Rlabel) + pf_type + ".txt");
  ifstream mc_file("./plots/indirectRho/" + pf_type + Form("/Spring16_25nsV1_MC_L1RC_AK%iPF", Rlabel) + pf_type + ".txt");
  string data_line, mc_line;

  //read first line
  getline(data_file, data_line);
  getline(mc_file, mc_line);

  const int ETA_BINS = 82;

  double eta1[ETA_BINS], eta2[ETA_BINS], data_p0[ETA_BINS], data_p1[ETA_BINS], data_p2[ETA_BINS];
  double mc_p0[ETA_BINS], mc_p1[ETA_BINS], mc_p2[ETA_BINS];

  for (int i=0; getline(data_file,data_line); i++){

    string str;
    int delim_pos;

    while (data_line.at(0) == ' ')  //positive values of eta have space in front of them
      data_line.erase(0, 1);

    //loop over strings in data line
    for (int col_num=0; (delim_pos = data_line.find(' ')) != -1; col_num++){

      str = data_line.substr(0, delim_pos);
      data_line.erase(0, delim_pos + 1);

      while (data_line.at(0) == ' ')  //get rid of white space between columns
        data_line.erase(0, 1);

      if (col_num == 0) eta1[i] = stod(str);
      else if (col_num == 1) eta2[i] = stod(str);
      else if (col_num == 9) data_p0[i] = stod(str);
      else if (col_num == 10) data_p1[i] = stod(str);
    }
    //last column after loop
    data_p2[i] = stod(data_line);

    //mc line
    getline(mc_file, mc_line);

    while (mc_line.at(0) == ' ')  //positive values of eta have space in front of them
      mc_line.erase(0, 1);

    //loop over strings in mc line
    for (int col_num=0; (delim_pos = mc_line.find(' ')) != -1; col_num++){

      str = mc_line.substr(0, delim_pos);
      mc_line.erase(0, delim_pos + 1);

      while (mc_line.at(0) == ' ')  //get rid of white space between columns
        mc_line.erase(0, 1);

      if (col_num == 9) mc_p0[i] = stod(str);
      else if (col_num == 10) mc_p1[i] = stod(str);
    }
    //last column after loop
    mc_p2[i] = stod(mc_line);
  }
  data_file.close();
  mc_file.close();

  TString hname;
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

  ofstream writeFile(Form("./plots/scalefactor/Spring16_25nsV1_DataMcSF_L1RC_AK%iPF", Rlabel) + pf_type + ".txt");
  writeFile << "{1   JetEta   1   Rho   [0]+[1]*x+[2]*pow(x,2)   Data/MC   L1FastJet}" << endl;

  TF1* fit = new TF1("fit", "1++x++x*x");
  fit->SetLineColor(1);
  fit->SetLineWidth(2);

  for (int i=0; i<ETA_BINS; i++){

    vector<double> rho, scale_factor, sf_error;

    for (int n=n1; n<n2; n++){

      double data_rho = data_rho_nPU->GetBinContent( data_rho_nPU->FindBin(n) );
      rho.push_back( data_rho );
      double data_offset = data_p0[i] + data_p1[i]*data_rho + data_p2[i]*data_rho*data_rho;

      double mc_rho = mc_rho_nPU->GetBinContent( mc_rho_nPU->FindBin(n) );
      double mc_offset = mc_p0[i] + mc_p1[i]*mc_rho + mc_p2[i]*mc_rho*mc_rho;

      scale_factor.push_back( data_offset / mc_offset );
      sf_error.push_back( 0.02*scale_factor.back() );
    }

    TGraphErrors* graph = new TGraphErrors(rho.size(), &rho[0], &scale_factor[0], 0, &sf_error[0]);
    graph->Fit(fit, "Q");

    writeFile << eta1[i] << setw(8) << eta2[i] << setw(8) << 5 << setw(6) << 0 << setw(6) << 200
              << setw(15) << fit->GetParameter(0) << setw(15) << fit->GetParameter(1) << setw(15) << fit->GetParameter(2) << endl;

    TCanvas* c = new TCanvas("c", "c", 600, 600);
    int start = rho[0]; int end = int(rho[rho.size()-1])+1;
    TH1F* h = new TH1F("h", "h", (end-start)*2, start, end);
    int topY = 2;

    h->GetXaxis()->SetTitle("#rho_{Data} (GeV)");
    h->GetYaxis()->SetTitle("Scale Factor");
    h->GetYaxis()->SetTitleOffset(1.05);
    h->GetYaxis()->SetRangeUser(0, topY);
    h->Draw();

    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(1);
    graph->Draw("psame");

    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);

    if (pf_type.EqualTo("all"))
      text.DrawLatex(0.17, 0.96, Form("AK%i PF %4.3f #leq #eta #leq %4.3f", Rlabel, eta1[i], eta2[i]) );
    else
      text.DrawLatex(0.17, 0.96, Form("AK%i PF%s %4.3f #leq #eta #leq %4.3f", Rlabel, pf_type.Data(), eta1[i], eta2[i]) );

    text.DrawLatex(0.2, 0.88, Form("#chi^{2}/ndof = %4.2f/%i", fit->GetChisquare(), fit->GetNDF() ) );
    text.DrawLatex(0.2, 0.84, Form("p0 = %4.3f #pm %4.3f", fit->GetParameter(0), fit->GetParError(0) ) );
    text.DrawLatex(0.2, 0.8, Form("p1 = %4.3f #pm %4.3f", fit->GetParameter(1), fit->GetParError(1) ) );
    text.DrawLatex(0.2, 0.76, Form("p2 = %4.4f #pm %4.4f", fit->GetParameter(2), fit->GetParError(2) ) );

    text.SetTextSize(0.035);
    text.SetTextColor(1);
    text.SetTextFont(42);
    text.DrawLatex(0.8, 0.96, "#sqrt{s} = 13 TeV");

    cout << fit->GetChisquare() / fit->GetNDF() << endl;

    c->Print("./plots/scalefactor/scalefactor_PF" + pf_type + Form("_eta%4.3f.pdf", eta1[i]) );
    delete h;
    delete c;
  }
  writeFile.close();
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
}
