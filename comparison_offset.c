// Bahareh Roozbahani and Chad Harrington Aug 2018
// This macro is for comparing offset pT derived from different methods, RC and FastJet
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include <map>

using namespace std;
void setStyle();
double cRC(const double& area, const double& rho, const double& pt, const double& p0, const double& p1, const double& p2) {
    double c = 1-(area/pt)*(p0+p1*(rho-1.519)+p2*pow(rho-1.519,2));
    return max(0.0001, c);
}
double cFastJet(const double& area, const double& rho, const double& pt, const double& p0, const double& p1, const double& p2) {
    double c = 1-area*(p0+(p1*rho)*(1+p2*log(pt)))/pt;
    return max(0.0001, c);
}
double cFastJet(const double& area, const double& rho, const double& pt, const double& p0, const double& p1, const double& p2, const double& p3, const double& p4, const double& p5) {
    double c = 1-(area/pt)*(p0+p1*(rho-20.0)+p2*log(pt/30.0)+p3*pow(log(pt/30.0),2)+p4*(rho-20.0)*log(pt/30.0)+p5*(rho-20.0)*pow(log(pt/30.0),2));
    return max(0.0001, c);
}

const int nEta = 82;
float etabins[nEta+1] =
{-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
    -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
    -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
    0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
    1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
    4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
void comparison_offset(int rlabel=4, TString pf_type="chs") {


    string run_name; float luminosity;
    cout << "Run: " << " lumi " << endl;
    cin >> run_name >> luminosity; 
   
    //Parameters
    bool versusPt = true; //versus pt or eta
    bool dependsOnMu = true;
    const double Rad=0.4;
    float eta = 0 , pt = 50;
    vector<float> mu = {30,40};
    vector<float> rho;
    
    TString plot_type;
    map<TString, TString> files;
    
    files["first"] = Form("Fall17_17Nov2017%s_V1_Data_L1FastJet_AK%iPF", run_name.c_str(), rlabel)+pf_type+".txt";
    //files["second"] = Form("Fall17_17Nov2017_V8_MC_L1FastJet_AK%iPF", rlabel)+pf_type+".txt";
    files["second"] = "plots/indirectRho/" + pf_type + Form("_%s/newFall17_17Nov2017%s_DATA_L1RC_AK%iPF", run_name.c_str(), run_name.c_str(), rlabel) + pf_type + ".txt";
    //files["second"] = "plots/indirectRho/" + pf_type + Form("_%s/newFall17_17Nov2017%s_DATA_L1RC_AK%iPF", run_name.c_str(), run_name.c_str(), rlabel) + pf_type + ".txt";
    
    map<TString, TString> rootfiles;
    map<TString, TString> files_label;
    rootfiles["first"] = Form("root_files_R48/Run%s_17nov_R%i.root",run_name.c_str(),rlabel);
    rootfiles["second"] = Form("root_files_R48/Run%s_17nov_R%i.root",run_name.c_str(),rlabel);
    
    
    TFile* rootfile = TFile::Open( rootfiles["second"] );
    TProfile* p_rho_nPU = (TProfile*) rootfile->Get("p_rho_nPU");
    for(int k=0 , n=mu.size(); k<n ; k++) rho.push_back( p_rho_nPU->GetBinContent( p_rho_nPU->FindBin( mu[k] ) ) );
    
    TH1F* h_nPU = (TH1F*) rootfile->Get("nPU");
    double mean_nPU = h_nPU->GetMean();
    //rho.push_back(p_rho_nPU->GetBinContent( p_rho_nPU->FindBin( mean_nPU ) ) );
    
    
    //End Parameters
    
    
    //plotting basics//
    setStyle();
    TCanvas* c = new TCanvas("c", "c", 600, 600);
    float b_scale = 0.3;
    float t_scale = 1 - b_scale;
    TPad* top = new TPad("top", "top", 0, b_scale, 1, 1);
    TPad* bottom = new TPad("bottom", "bottom", 0, 0, 1, b_scale);
    top->SetTopMargin(0.05);
    top->SetBottomMargin(0.05);
    top->Draw();
    bottom->SetBottomMargin(0.35);
    bottom->SetGridy();
    bottom->Draw();
    top->cd();
    
    
    TH1F* h;
    int bins = 100, xmax = 1000;
    if (versusPt) {
        h = new TH1F("h", "h", bins, 0, xmax);
        h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
        h->GetXaxis()->SetNdivisions(5, 5, 0);
    }
    else {
        h = new TH1F("h", "h", nEta, etabins);
        h->GetXaxis()->SetTitle("#eta");
    }
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetTitle("Offset p_{T}");
    h->Draw();
    TLegend* leg = new TLegend(.7,.8-(files.size()*0.05),.9,.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    
    TLegend* leg2 = new TLegend(.2,.8-(files.size()*0.05),.5,.85);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
    leg2->SetTextFont(42);
    // plotting basics//
    
    int etabin = -1;
    if (versusPt) {
        
        for (int i=0; i<nEta+1; i++) {
            if (etabins[i] <= eta && eta <etabins[i+1]) { etabin = i; break; }
        }
        if (etabin == -1) {cout << "eta out of range" << endl; return; }
    }
    
    
    map<TString, double[nEta]> p0, p1, p2, p3, p4, p5;
    map< pair<TString, float> , TH1F*> hists;
    
    for (map<TString, TString>::iterator it = files.begin(); it != files.end(); it++){  
      TString key = it->first;
      cout<<key<<endl;
      ifstream file(it->second);
      string line;
        
      //read first line
      getline(file, line);
      TString mc_formula = TString(line) ;
      TString mapkey = "";
      if (!mc_formula.Contains("pow(log(y/30.0),2") && !mc_formula.Contains("log(y/90.0)") ){
        cout << "old methods" << endl;    
        cout << line << endl;
        if (mc_formula.Contains("1-z*([0]+([1]*x)*(1+[2]*log(y)))/y") || mc_formula.Contains("1-y*([1]*(z-[0])*(1+[2]*log(x/30.)))/x") ) {
          mapkey = "oldFastjet";
        }
        else { mapkey = "L1RC";}

        for (int i=0; getline(file,line); i++){    
          string str;
          int delim_pos;    
          while (line.at(0) == ' ')  //eta values have space in front of them
            line.erase(0, 1);
            
          //loop over strings in line
          for (int col_num=0; (delim_pos = line.find(' ')) != -1; col_num++){        
            str = line.substr(0, delim_pos);
            line.erase(0, delim_pos + 1);
                
            while (line.at(0) == ' ')  //get rid of white space between columns
              line.erase(0, 1);
                
            if (col_num == 9) p0[key][i] = stod(str);
            else if (col_num == 10) p1[key][i] = stod(str);
          }
          //last column after loop
          p2[key][i] = stod(line);
          //cout << mapkey << "\t" << p0[key][i] << "\t" << p1[key][i] << "\t" << p2[key][i] << endl;
        }
      }
      else{
        cout << "new method" << endl;
        int lines;
        file.clear();
        file.seekg(0, file.beg);
        getline(file,line);
        cout << line <<endl;
        mapkey = "newFastjet";
        for (int i=0; getline(file,line); i++){
          lines = i;
          string str;
          int delim_pos;
          while (line.at(0) == ' ')  //eta values have space in front of them
          line.erase(0, 1);
          //loop over strings in mc line
          delim_pos = line.find(' ');
          int it = -10;
          for (int col_num=0; (delim_pos = line.find(' ')) != -1 ; col_num++){
            str = line.substr(0, delim_pos);
            line.erase(0, delim_pos + 1);

            while (!line.empty() && line.at(0) == ' ')  //get rid of white space between columns
              line.erase(0, 1);

            if (col_num == 9) p0[key][i] = stod(str);
            else if (col_num == 10) p1[key][i] = stod(str);
            else if (col_num == 11) p2[key][i] = stod(str);
            else if (col_num == 12) p3[key][i] = stod(str);
            else if (col_num == 13) p4[key][i] = stod(str);
            it = col_num;
          }
          if ((delim_pos = line.find(' ')) == -1 && it == 13) p5[key][i] = stod(line.substr(0, line.length()));
          //cout << mapkey << "\t" << p0[key][i] << "\t" << p1[key][i] << "\t" <<  p2[key][i] << "\t" << p3[key][i] << "\t" << p4[key][i] << "\t" << p5[key][i] << "\t" << endl;
        }
      } 
      file.close();        
      int j = distance( files.begin(), it );
      float area = M_PI*Rad*Rad;
      double x = fabs(etabins[0]) - 0.5*fabs(etabins[etabin]+etabins[etabin]);
      if (x<Rad){
        double theta = 2*acos(x/Rad);
        double area_seg = 0.5*Rad*Rad*theta - x*Rad*sin(theta/2);
         area -= area_seg;
      }
      if (versusPt) {
        plot_type = "_pT";
        for(int k=0 , n=rho.size(); k<n ; k++) hists[make_pair(key,rho[k])] = new TH1F( Form("h%i%i",j,k), Form("h%i%i",j,k), bins, 0, xmax );
        for(int k=0 , n=rho.size(); k<n ; k++){
          for (int i=0; i<bins; i++) {
            float i_pt = (i+1)*xmax/bins;
            double c;
            if (mapkey.Contains("L1RC", TString::kIgnoreCase)){
              c = cRC(area, rho[k], i_pt, p0[key][etabin], p1[key][etabin], p2[key][etabin]);

            }
            else if (mapkey.Contains("Fastjet", TString::kIgnoreCase)){
              if (mapkey.Contains("old", TString::kIgnoreCase)) c = cFastJet(area, rho[k], i_pt, p0[key][etabin], p1[key][etabin], p2[key][etabin]);
              else if (mapkey.Contains("new", TString::kIgnoreCase)) c = cFastJet(area, rho[k], i_pt, p0[key][etabin], p1[key][etabin], p2[key][etabin]
, p3[key][etabin], p4[key][etabin], p5[key][etabin]);
            }
            hists[make_pair(key,rho[k])]->SetBinContent(i+1, i_pt-c*i_pt);
            cout << mapkey << "\t" << i_pt << "\t" << c << "\t" << i_pt-c*i_pt << endl;
          }
        }
      }
      else {
        plot_type = "_eta";
        for(int k=0 , n=rho.size(); k<n ; k++) hists[make_pair(key,rho[k])] = new TH1F( Form("h%i%i",j,k), Form("h%i%i",j,k), nEta, etabins);
        for(int k=0 , n=rho.size(); k<n ; k++){
          for (int i_eta=0; i_eta<nEta; i_eta++) {
            double c;
            float area = M_PI*Rad*Rad;
            double x = fabs(etabins[0]) - 0.5*fabs(etabins[i_eta]+etabins[i_eta+1]);
            //if (key.Contains("RC", TString::kIgnoreCase))
            if (x<Rad){
              double theta = 2*acos(x/Rad);
              double area_seg = 0.5*Rad*Rad*theta - x*Rad*sin(theta/2);
              area -= area_seg;
            }
            if (mapkey.Contains("L1RC", TString::kIgnoreCase)){
              c = cRC(area, rho[k], pt, p0[key][i_eta], p1[key][i_eta], p2[key][i_eta]);
            }
            else if (mapkey.Contains("Fastjet", TString::kIgnoreCase)){
              if (mapkey.Contains("old", TString::kIgnoreCase)) c = cFastJet(area, rho[k], pt, p0[key][i_eta], p1[key][i_eta], p2[key][i_eta]);
              else if (mapkey.Contains("new", TString::kIgnoreCase)) c = cFastJet(area, rho[k], pt, p0[key][i_eta], p1[key][i_eta], p2[key][i_eta]
, p3[key][i_eta], p4[key][i_eta], p5[key][i_eta]);
            }
            hists[make_pair(key,rho[k])]->SetBinContent(i_eta+1, pt-c*pt);
          }
        }
      }

        
      hists[make_pair(key,rho[0])]->SetMarkerColor(1);
      hists[make_pair(key,rho[1])]->SetMarkerColor(2);
      //hists[make_pair(key,rho[2])]->SetMarkerColor(4);
      //hists[make_pair(key,rho[3])]->SetMarkerColor(8);
      hists[make_pair(key,rho[0])]->SetLineColor(1);
      hists[make_pair(key,rho[1])]->SetLineColor(2);
      //hists[make_pair(key,rho[2])]->SetLineColor(4);
      //hists[make_pair(key,rho[3])]->SetLineColor(8);
      for(int k=0 , n=rho.size(); k<n ; k++){
        hists[make_pair(key,rho[k])]->SetMarkerSize(0.75);
        hists[make_pair(key,rho[k])]->Draw("psame");        
      }
      //leg2->AddEntry(hists[make_pair(key,15)], key ,"L");
      if (mapkey.Contains("L1RC")) files_label[key] = "RC";
      else if (mapkey.Contains("OLDFast",TString::kIgnoreCase)) files_label[key] = "V8_FastJet";
      else if (mapkey.Contains("NEWFast",TString::kIgnoreCase)) files_label[key] = "V22_FastJet";
    }
    for(int i=0, n=rho.size(); i<n ; i++){
      hists[make_pair("first",rho[i])]->SetMarkerStyle(4);
      hists[make_pair("second",rho[i])]->SetMarkerStyle(22);
    }
    leg2->AddEntry(hists[make_pair("first",rho[0])], files_label["first"] ,"P");
    leg2->AddEntry(hists[make_pair("second",rho[0])], files_label["second"] ,"P");
    for(int i=0, n=rho.size(); i<n; i++) leg->AddEntry(hists[make_pair("first",rho[i])], Form("#rho = %2.0f",rho[i]),"L");
    //leg->AddEntry(hists[make_pair("first",rho[3])], Form("nominal #rho = %2.0f",rho[3]),"L");
    leg->Draw();
    leg2->Draw();
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.04);
    
    if (versusPt)
        text.DrawLatex( 0.25, 0.87, Form("%1.3f#leq#eta<%1.3f", etabins[etabin], etabins[etabin+1]) );
    else
        text.DrawLatex( 0.25, 0.87, Form("p_{T} = %2.0f GeV", pt) );
        
        if (pf_type.EqualTo("all"))
            text.DrawLatex(0.17, 0.96, Form("AK%i PF", rlabel));
        else
            text.DrawLatex(0.17, 0.96, Form("AK%i PF%s", rlabel, pf_type.Data()));
    
    text.SetTextSize(0.035);
    text.SetTextColor(1);
    text.SetTextFont(42);
    text.DrawLatex(0.72, 0.96, Form("2017%s (13 TeV)",run_name.c_str()));
    //text.DrawLatex(0.72, 0.96, Form("MC",run_name.c_str()));
    if (files_label["first"].Contains("v22",TString::kIgnoreCase) && files_label["second"].Contains("RC",TString::kIgnoreCase)) h->GetYaxis()->SetRangeUser(-6.2,10.2);
    if (files_label["first"].Contains("v8",TString::kIgnoreCase) && files_label["second"].Contains("RC",TString::kIgnoreCase)) h->GetYaxis()->SetRangeUser(3,8.2);
    if (files_label["first"].Contains("v22",TString::kIgnoreCase) && files_label["second"].Contains("V8",TString::kIgnoreCase)) h->GetYaxis()->SetRangeUser(-3,16.2);
    
    bottom->cd();
    TH1D* h2;
    if (versusPt) h2 = new TH1D("h2", "h2", bins, 0, xmax);
    else h2 = new TH1D("h2", "h2", 82, etabins);
    
    map<pair<TString, float> ,TH1D*> first;
    map<pair<TString, float> ,TH1D*> second;
    for(int i=0, n=rho.size();i<n;i++){
       first[make_pair("first",rho[i])] = (TH1D*) hists[make_pair("first",rho[i])]->Clone(Form("first_%f",rho[i]));
       second[make_pair("second",rho[i])] = (TH1D*) hists[make_pair("second",rho[i])]->Clone(Form("second_%f",rho[i]));
    }
    for(int i=0, n=rho.size();i<n;i++) first[make_pair("first",rho[i])]->Divide(second[make_pair("second",rho[i])]);
    
    
    h2->GetXaxis()->SetLabelSize(0.03/b_scale);
    h2->GetXaxis()->SetTickLength(0.03/b_scale);
    h2->GetXaxis()->SetTitleSize(0.035/b_scale);
    h2->GetXaxis()->SetTitleOffset(1);
    if (versusPt) {
        h2->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
    }
    else h2->GetXaxis()->SetTitle("#eta");
    h2->GetYaxis()->SetNdivisions(5, 3, 0);
    h2->GetYaxis()->SetLabelSize(0.025/b_scale);
    h2->GetYaxis()->SetTitle(files_label["first"]+"/"+files_label["second"]);
    h2->GetYaxis()->SetTitleSize(0.035/b_scale);
    h2->GetYaxis()->SetTitleOffset(0.43);
    
    for(int i=0, n=rho.size();i<n;i++)  first[make_pair("first",rho[i])]->SetMarkerStyle(4);
    first[make_pair("first",rho[0])]->SetMarkerColor(1);
    first[make_pair("first",rho[1])]->SetMarkerColor(2);
    //first[make_pair("first",rho[2])]->SetMarkerColor(4);
    //first[make_pair("first",rho[3])]->SetMarkerColor(8);
    h2->Draw();
    for(int i=0, n=rho.size();i<n;i++)  first[make_pair("first",rho[i])]->Draw("sameP");
    if (files_label["first"].Contains("v22",TString::kIgnoreCase) && files_label["second"].Contains("RC",TString::kIgnoreCase)) h2->GetYaxis()->SetRangeUser(-2,2);
    if (files_label["first"].Contains("v8",TString::kIgnoreCase) && files_label["second"].Contains("RC",TString::kIgnoreCase)) h2->GetYaxis()->SetRangeUser(0.9,1.3);
    TLegend* leg3 = new TLegend(.7,.8-(files.size()*0.05),.9,.9);
    leg3->SetBorderSize(0);
    leg3->SetFillColor(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.06);
    leg3->SetTextFont(42);
    //for(int i=0, n=rho.size()-1;i<n;i++)  leg3->AddEntry(SF[make_pair("Data_L1RC_SF",rho[i])],Form("#rho = %1.0f",rho[i]),"P");
    //leg3->AddEntry(first[make_pair("first",rho[3])],Form("nominal #rho = %1.0f",rho[3]),"P");
    leg3->Draw();
    
    
    
    c->Print("offset_pt_PF" + pf_type + plot_type + ".pdf");
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
