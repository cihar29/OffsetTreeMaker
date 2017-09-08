//Chad Harrington 6/14/2017
//EXECUTE as root -l -b -q 'offsetpT_stack.c ("MC_R4.root", "Data_R4.root", "nPU", Id, ratio, "label")'

using namespace std;

void setStyle();

const int ETA_BINS = 82;
float etabins[ETA_BINS+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

enum Id{
  ne=0, hfe, nh, hfh, chu, chm, untrk, numId, all, hf_dep
};

void offsetpT_stack( TString mcName="histoMC_R4.root", TString dataName="histoData_R4.root", TString bin_var="nPU", Id id=all, bool ratio=true,
                     TString label="Run 2016 - 35.9 fb^{-1} (13 TeV)") {

  TFile* mcFile = TFile::Open(mcName);
  TFile* dataFile = TFile::Open(dataName);

  TString ids[] = {"ne", "hfe", "nh", "hfh", "chu", "chm", "untrk", "NumIds", "all", "hf_dep"};

  //nPV or nPU
  TH1F* h_bin_var = (TH1F*) dataFile->Get(bin_var);
  //int n1 = h_bin_var->GetBinCenter( h_bin_var->GetMaximumBin() );
  int n1 = h_bin_var->GetMean();

  vector<TH1D*> v_MC (numId);
  vector<TH1D*> v_Data (numId);

  TString hname;
  for (int i=0; i<numId; i++) {

    hname = Form("p_offset_eta_%s%i_", bin_var.Data(), n1) + ids[i];

    v_MC[i] = ((TProfile*) mcFile->FindObjectAny(hname))->ProjectionX(ids[i]+"MC");
    v_Data[i] = ((TProfile*) dataFile->FindObjectAny(hname))->ProjectionX(ids[i]+"Data");

    ///Regular Scaling///
    v_MC[i]->Scale( 1.0 / n1 );
    v_Data[i]->Scale( 1.0 / n1 );
  }

  /////////////////////////////////////////
  /// Slightly More Complicated Scaling ///
  /////////////////////////////////////////
/*
  int n2 = 50;
  for (int i=0; i<numIds; i++) {

    v_MC[i]->Scale(-1);
    v_Data[i]->Scale(-1);

    hname = Form("p_offset_eta_%s%i_", bin_var.Data(), n2) + ids[i];

    v_MC[i]->Add( ((TProfile*) mcFile->FindObjectAny(hname))->ProjectionX(ids[i]+"MC2") );
    v_Data[i]->Add( ((TProfile*) dataFile->FindObjectAny(hname))->ProjectionX(ids[i]+"Data2") );

    v_MC[i]->Scale( 1.0 / (n2-n1) );
    v_Data[i]->Scale( 1.0 / (n2-n1) );
  }
*/
  ///////////////////////////
  /// Complicated Scaling ///
  ///////////////////////////
/*
  TH1D* h_weight_MC = (TH1D*) mcFile->Get(bin_var);
  TH1D* h_weight_Data = (TH1D*) dataFile->Get(bin_var);

  if (h_weight_Data->GetXaxis()->GetBinWidth(1) == 0.5) {
    h_weight_MC->Rebin();
    h_weight_Data->Rebin();
  }
  h_weight_MC->Scale( 1 / h_weight_MC->Integral() );
  h_weight_Data->Scale( 1 / h_weight_Data->Integral() );

  for (int i=0; i<numIds; i++) {

    double total_weight_MC = 0;
    double total_weight_Data = 0;

    TH1D* sum_MC = new TH1D("sum_MC", "sum_MC", ETA_BINS, etabins);
    TH1D* sum_Data = new TH1D("sum_Data", "sum_Data", ETA_BINS, etabins);

    for (int j=n1+1; j<=n2; j++) {

      hname = Form("p_offset_eta_%s%i_", bin_var.Data(), j) + ids[i];

      TH1D* temp_MC = ((TProfile*) mcFile->FindObjectAny(hname))->ProjectionX("temp_MC");
      TH1D* temp_Data = ((TProfile*) dataFile->FindObjectAny(hname))->ProjectionX("temp_Data");

      temp_MC->Add(v_MC[i], -1);
      temp_Data->Add(v_Data[i], -1);

      double weight_MC = h_weight_MC->GetBinContent( h_weight_MC->FindBin(j) );
      double weight_Data = h_weight_Data->GetBinContent( h_weight_Data->FindBin(j) );

      temp_MC->Scale( weight_MC / (j-n1) );
      temp_Data->Scale( weight_Data / (j-n1) );

      sum_MC->Add(temp_MC);
      sum_Data->Add(temp_Data);

      delete temp_MC;
      delete temp_Data;

      total_weight_MC += weight_MC;
      total_weight_Data += weight_Data;
    }
    v_MC[i] = (TH1D*) sum_MC->Clone(ids[i]);
    v_Data[i] = (TH1D*) sum_Data->Clone(ids[i]);

    delete sum_MC;
    delete sum_Data;

    v_MC[i]->Scale( 1 / total_weight_MC );
    v_Data[i]->Scale( 1 /total_weight_Data );
  }
*/
  //Make PF all and PF chs for ratio plot before alterations//

  TH1D* all_Data = (TH1D*) v_Data[untrk-1]->Clone("all_Data");  //don't include lost tracks
  TH1D* all_MC = (TH1D*) v_MC[untrk-1]->Clone("all_MC");

  // ne + hfe + nh + hfh + chu + chm
  for (int i=0; i<untrk-1; i++){
    all_Data->Add(v_Data[i]);
    all_MC->Add(v_MC[i]);
  }

  //chs is everything except matched hadrons
  TH1D* chs_Data = (TH1D*) v_Data[chu]->Clone("chs_Data");
  TH1D* chs_MC = (TH1D*) v_MC[chu]->Clone("chs_MC");

  //ne + hfe + nh + hfh + chu
  for (int i=0; i<chu; i++){
    chs_Data->Add(v_Data[i]);
    chs_MC->Add(v_MC[i]);
  }

  //Draw Markers for EM Deposits and Hadronic Deposits in two separate regions//

  TH1D* EM_clone = (TH1D*) v_Data[hfe]->Clone("EM_clone");
  TH1D* Had_clone = (TH1D*) v_Data[hfh]->Clone("Had_clone");

  //Stacks//

  THStack* mcStack = new THStack();
  THStack* dataStack = new THStack();
  THStack* cloneStack = new THStack();

  for (int i=0; i<numId-1; i++){  //don't add lost tracks
    mcStack->Add(v_MC[i]);
    dataStack->Add(v_Data[i]);
  }

  cloneStack->Add(v_Data[ne]);
  cloneStack->Add(EM_clone);
  cloneStack->Add(v_Data[nh]);
  cloneStack->Add(Had_clone);

  setStyle();

  TString yTitle = "";
  if ( bin_var.EqualTo("nPV") ) yTitle = "<Offset p_{T}> / <N_{PV}> (GeV)";
  else yTitle = "<Offset p_{T}> / <#mu> (GeV)";

  TCanvas* c = new TCanvas("c", "c", 600, 600);
  gStyle->SetOptStat(0);

  float b_scale = 0.3, t_scale = 1 - b_scale;

  TPad* top = new TPad("top", "top", 0, b_scale, 1, 1);
  TPad* bottom = new TPad("bottom", "bottom", 0, 0, 1, b_scale);

  if (ratio) {
    top->SetTopMargin(0.05);
    top->SetBottomMargin(0.05);
    top->Draw();
    bottom->SetBottomMargin(0.35);
    bottom->SetGridy();
    bottom->Draw();
    top->cd();
  }

  TH1D* h1 = new TH1D("h1", "h1", ETA_BINS, etabins);
  h1->GetYaxis()->SetTitle(yTitle);

  if (ratio) {
    h1->GetXaxis()->SetTickLength(0.03/t_scale);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetYaxis()->SetTitleSize(0.05/t_scale);
    h1->GetYaxis()->SetTitleOffset(0.9);
    h1->GetYaxis()->SetLabelSize(0.04/t_scale);
  }
  else {
    h1->GetXaxis()->SetTitle("#eta");
    h1->GetXaxis()->SetTitleSize(0.03/t_scale);
    h1->GetXaxis()->SetTitleOffset(1.1);
    h1->GetXaxis()->SetLabelSize(0.03/t_scale);

    h1->GetYaxis()->SetTitleSize(0.03/t_scale);
    h1->GetYaxis()->SetTitleOffset(1.3);
    h1->GetYaxis()->SetLabelSize(0.03/t_scale);
  }
  h1->Draw();

  v_MC[ne]   ->SetMarkerStyle(kMultiply);
  v_MC[hfe]  ->SetMarkerStyle(kOpenStar);
  v_MC[nh]   ->SetMarkerStyle(kOpenDiamond);
  v_MC[hfh]  ->SetMarkerStyle(kOpenTriangleUp);
  v_MC[chu]  ->SetMarkerStyle(kOpenCircle);
  v_MC[chm]  ->SetMarkerStyle(kOpenCircle);
  v_MC[untrk]->SetMarkerStyle(kOpenCircle);

  v_Data[ne]    ->SetMarkerStyle(kMultiply);
  v_Data[hfe]   ->SetMarkerStyle(kOpenStar);
  EM_clone      ->SetMarkerStyle(kOpenStar);
  v_Data[nh]    ->SetMarkerStyle(kOpenDiamond);
  v_Data[hfh]   ->SetMarkerStyle(kOpenTriangleUp);
  Had_clone     ->SetMarkerStyle(kOpenTriangleUp);
  v_Data[chu]   ->SetMarkerStyle(kOpenCircle);
  v_Data[chm]   ->SetMarkerStyle(kOpenCircle);
  v_Data[untrk] ->SetMarkerStyle(kOpenCircle);

  //for error bars
  v_Data[ne]    ->SetLineColor(kBlack);
  v_Data[hfe]   ->SetLineColor(kBlack);
  EM_clone      ->SetLineColor(kBlack);
  v_Data[nh]    ->SetLineColor(kBlack);
  v_Data[hfh]   ->SetLineColor(kBlack);
  Had_clone     ->SetLineColor(kBlack);
  v_Data[chu]   ->SetLineColor(kBlack);
  v_Data[chm]   ->SetLineColor(kBlack);
  v_Data[untrk] ->SetLineColor(kBlack);

  v_MC[ne]   ->SetFillColor(kBlue);
  v_MC[hfe]  ->SetFillColor(kViolet+2);
  v_MC[nh]   ->SetFillColor(kGreen);
  v_MC[hfh]  ->SetFillColor(kPink+6);
  v_MC[chu]  ->SetFillColor(kRed-9);
  v_MC[chm]  ->SetFillColor(kRed);
  v_MC[untrk]->SetFillColor(kOrange);

  v_MC[ne]   ->SetLineColor(kBlack);
  v_MC[hfe]  ->SetLineColor(kBlack);
  v_MC[nh]   ->SetLineColor(kBlack);
  v_MC[hfh]  ->SetLineColor(kBlack);
  v_MC[chu]  ->SetLineColor(kBlack);
  v_MC[chm]  ->SetLineColor(kBlack);
  v_MC[untrk]->SetLineColor(kBlack);

  //CMS and Lumi Text//

  TLatex text;
  text.SetNDC();

  if (id == all) {

    v_Data[ne]   ->SetAxisRange(-2.9,2.9);
    v_Data[hfe]  ->SetAxisRange(-5,-2.6);
    EM_clone     ->SetAxisRange(2.6,5);
    v_Data[hfh]  ->SetAxisRange(-5,-2.6);
    Had_clone    ->SetAxisRange(2.6,5);
    v_Data[nh]   ->SetAxisRange(-2.9,2.9);
    v_Data[chu]  ->SetAxisRange(-2.9,2.9);
    v_Data[chm]  ->SetAxisRange(-2.9,2.9);
    v_Data[untrk]->SetAxisRange(-2.9,2.9);

    h1->GetYaxis()->SetRangeUser( 0, 0.8 ); //dataStack->GetMaximum()*1.7 );
    mcStack->Draw("samehist");
    dataStack->Draw("samepe");
    cloneStack->Draw("samepe");

    TLegend* leg = new TLegend(.4,.57,.67,.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetHeader("#bf{Markers: Data, Histograms: MC}");
    leg->AddEntry(v_MC[ne],"Photons","PF");
    leg->AddEntry(v_MC[hfe],"EM Deposits","PF");
    leg->AddEntry(v_MC[nh],"Neutral Hadrons","PF");
    leg->AddEntry(v_MC[hfh],"Hadronic Deposits","PF");
    leg->AddEntry(v_MC[chu],"Unassoc. Charged Hadrons","PF");
    leg->AddEntry(v_MC[chm],"Assoc. Charged Hadrons","PF");
    //leg->AddEntry(v_MC[untrk],"Lost Tracks","PF");
    leg->Draw();

    text.SetTextSize(0.065);
    text.SetTextFont(61);
    text.DrawLatex(0.22, 0.85, "CMS");

    text.SetTextSize(0.045);
    text.SetTextFont(42);

    if (ratio) text.DrawLatex(1-label.Length()/80., 0.96, label);
    else       text.DrawLatex(1-label.Length()/68., 0.96, label);

    //TString coneSize = dataName( dataName.Last('.')-1, 1 );
    //text.DrawLatex(0.2, 0.8, "R = 0." + coneSize);

    gPad->RedrawAxis();

    if (ratio) {
      bottom->cd();
      TH1D* h2 = new TH1D("h2", "h2", ETA_BINS, etabins);

      chs_Data->Divide(chs_MC);
      all_Data->Divide(all_MC);

      h2->GetXaxis()->SetLabelSize(0.04/b_scale);
      h2->GetXaxis()->SetTickLength(0.03/b_scale);
      h2->GetXaxis()->SetTitleSize(0.05/b_scale);
      h2->GetXaxis()->SetTitleOffset(0.8);
      h2->GetXaxis()->SetTitle("#eta");
      h2->GetYaxis()->SetRangeUser(0.5, 2); //chs_Data->GetMaximum()*1.1 );
      h2->GetYaxis()->SetNdivisions(5, 3, 0);
      h2->GetYaxis()->SetLabelSize(0.04/b_scale);
      h2->GetYaxis()->SetTitle("Data/MC");
      h2->GetYaxis()->CenterTitle(true);
      h2->GetYaxis()->SetTitleSize(0.05/b_scale);
      h2->GetYaxis()->SetTitleOffset(0.4);

      chs_Data->SetLineColor(kBlack);
      all_Data->SetLineColor(kBlack);

      chs_Data->SetMarkerStyle(24);
      chs_Data->SetMarkerColor(2);
      all_Data->SetMarkerStyle(24);
      h2->Draw();
      chs_Data->Draw("samePE");
      all_Data->Draw("samePE");

      TLegend* leg2 = new TLegend(.55,.75,.65,.9);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetFillStyle(0);
      leg2->SetTextSize(0.06);
      leg2->SetTextFont(42);
      leg2->AddEntry(chs_Data,"PF chs","P");
      leg2->AddEntry(all_Data,"PF","P");
      leg2->Draw();
    }
  }
  else {
    TH1D* hist_MC;
    TH1D* hist_Data;

    if (id != hf_dep){
      hist_MC = v_MC[id];
      hist_Data = v_Data[id];
    }
    else{
      //HF = EM Deposits + Hadronic Deposits
      hist_MC = (TH1D*) v_MC[hfe]->Clone("hist_MC");
      hist_MC->Add( v_MC[hfh] );
      hist_Data = (TH1D*) v_Data[hfe]->Clone("hist_Data");
      hist_Data->Add( v_Data[hfh] );
    }

    h1->GetYaxis()->SetRangeUser( 0, 0.5 ); //hist_Data->GetMaximum()*1.4 );

    TString title;
    if (id == ne) { title = "Photons"; hist_Data->SetAxisRange(-2.9, 2.9); }
    else if (id == hfe) title = "EM Deposits";
    else if (id == nh) { title = "Neutral Hadrons"; hist_Data->SetAxisRange(-2.9, 2.9); }
    else if (id == hfh) title = "Hadronic Deposits";
    else if (id == chu) { title = "Unassoc. Charged Hadrons"; hist_Data->SetAxisRange(-2.9, 2.9); }
    else if (id == chm) { title = "Assoc. Charged Hadrons"; hist_Data->SetAxisRange(-2.9, 2.9); }
    else if (id == untrk) { title = "Lost Tracks"; hist_Data->SetAxisRange(-2.9, 2.9); }
    else title = "HF Deposits";

    hist_MC->Draw("sameHIST");
    hist_Data->Draw("samePE");

    TLegend* leg = new TLegend(.7,.7,.9,.8);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->AddEntry(hist_Data,"Data","P");
    leg->AddEntry(hist_MC,"MC","F");
    leg->Draw();

    text.SetTextSize(0.035/t_scale);
    text.SetTextFont(61);
    text.DrawLatex(0.17, 0.96, title);

    text.SetTextSize(0.065);
    text.SetTextFont(61);
    text.DrawLatex(0.25, 0.8, "CMS");

    text.SetTextSize(0.045);
    text.SetTextFont(42);

    if (ratio) text.DrawLatex(1-label.Length()/80., 0.96, label);
    else       text.DrawLatex(1-label.Length()/68., 0.96, label);

    //TString coneSize = dataName( dataName.Last('.')-1, 1 );
    //text.DrawLatex( 0.25, 0.75, "R = 0." + coneSize );

    gPad->RedrawAxis();

    if (ratio){

      bottom->cd();
      TH1D* h2 = new TH1D("h2", "h2", ETA_BINS, etabins);

      TH1D* ratio_MC = (TH1D*) hist_MC->Clone("ratio_MC");
      TH1D* ratio_Data = (TH1D*) hist_Data->Clone("ratio_Data");
      ratio_Data->Divide(ratio_MC);

      if (id == ne) ratio_Data->SetAxisRange(-2.9, 2.9);
      //else if (id == hfe) //em deposits
      else if (id == nh) ratio_Data->SetAxisRange(-2.9, 2.9);
      //else if (id == hfh) //hadronic deposits
      else if (id == chu) ratio_Data->SetAxisRange(-2.9, 2.9);
      else if (id == chm) ratio_Data->SetAxisRange(-2.9, 2.9);
      else if (id == untrk) ratio_Data->SetAxisRange(-2.9, 2.9);
      else if (id == hf_dep) ratio_Data->SetAxisRange(-3, 3);
      //else

      h2->GetXaxis()->SetLabelSize(0.04/b_scale);
      h2->GetXaxis()->SetTickLength(0.03/b_scale);
      h2->GetXaxis()->SetTitleSize(0.05/b_scale);
      h2->GetXaxis()->SetTitleOffset(0.8);
      h2->GetXaxis()->SetTitle("#eta");
      h2->GetYaxis()->SetRangeUser(0, 2);
      h2->GetYaxis()->SetNdivisions(5, 3, 0);
      h2->GetYaxis()->SetLabelSize(0.04/b_scale);
      h2->GetYaxis()->SetTitle("Data/MC");
      h2->GetYaxis()->CenterTitle(true);
      h2->GetYaxis()->SetTitleSize(0.05/b_scale);
      h2->GetYaxis()->SetTitleOffset(0.4);

      ratio_Data->SetMarkerStyle(24);
      h2->Draw();
      ratio_Data->Draw("samePE");
    }
  }

  c->Print("plots/stack_" + ids[id] + "_" + bin_var + to_string(n1) + ".pdf");
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

  tdrStyle->cd();
}
