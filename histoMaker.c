//Chad Harrington 6/8/2015
//EXECUTE as root -l -b -q histoMaker.c

#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <cmath>
#include <iomanip>
#include <vector>
using namespace std;

const int nEta = 82;
float etabins[nEta+1] =
  {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65,
   -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
   -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0,
   0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479,
   1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013,
   4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

void FillHist1D(const TString& histName, const Double_t& value, double weight);
void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, double weight);
void FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, double weight);
void FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, double weight);
void getGeometry(double (&geo)[nEta][nEta]);
double areaS(double R, double x1, double x2);
double dist(double R, double x1, double x2);

map<TString, TH1*> m_Histos1D;
map<TString, TH2*> m_Histos2D;
map<TString, TProfile*> m_Profiles;
map<TString, TProfile2D*> m_Profiles2D;

const int MAXNPU = 50;
const int MAXNPV = 50;
const int MAXRHO = 50;
const bool ISMC = false;
const bool REWEIGHT = true;
const float rCone = 0.4;

void histoMaker(){

  //Open Files//

  TString OUTFILENAME = Form("Data76x_R%i.root", int(rCone*10));
  TFile* inFile = TFile::Open("Data101_76x.root");

  TTree* T = (TTree*) inFile->Get("T");
  Long64_t nEntries = T->GetEntries();
  cout << nEntries << " Events" << endl;

  //Declare Histos//

  enum Flavor{ chm = 0, chu, nh, ne, hfh, hfe, lep, untrk, numFlavors};
  TString ids[] = {"chm", "chu", "nh", "ne", "hfh", "hfe", "lep", "untrk"};
  TString hname;

  for (int i_id=0; i_id<numFlavors; i_id++){
    for (int i_nPU=0; i_nPU<MAXNPU; i_nPU++){
      hname = Form("p_offset_eta_nPU%i_", i_nPU) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, nEta, etabins);
    }
    for (int i_nPV=0; i_nPV<MAXNPV; i_nPV++){
      hname = Form("p_offset_eta_nPV%i_", i_nPV) + ids[i_id];
      m_Profiles[hname] = new TProfile(hname, hname, nEta, etabins);
    }
  }

  //hname = "nPV_all";
  //m_Histos1D[hname] = new TH1F(hname,hname,MAXNPV,0,MAXNPV);
  hname = "nPV";
  m_Histos1D[hname] = new TH1F(hname,hname,MAXNPV,0,MAXNPV);
  //hname = "nPV_ratio";
  //m_Histos1D[hname] = new TH1F(hname,hname,100,0,1);
  //hname = "pv_z";
  //m_Histos1D[hname] = new TH1F(hname,hname,100,-50,50);
  hname = "nPU";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,MAXNPU);
  hname = "rho";
  m_Histos1D[hname] = new TH1F(hname,hname,100,0,MAXRHO);

  //hname = "2h_nPV_nPU";
  //m_Histos2D[hname] = new TH2F(hname,hname,100,0,MAXNPU,MAXNPV,0,MAXNPV);
  hname = "p_nPV_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,MAXNPU);

  //hname = "2h_rho_nPU";
  //m_Histos2D[hname] = new TH2F(hname,hname,100,0,MAXNPU,100,0,MAXRHO);
  hname = "p_rho_nPU";
  m_Profiles[hname] = new TProfile(hname,hname,100,0,MAXNPU);

  //hname = "2h_rho_nPV";
  //m_Histos2D[hname] = new TH2F(hname,hname,MAXNPV,0,MAXNPV,100,0,MAXRHO);
  hname = "p_rho_nPV";
  m_Profiles[hname] = new TProfile(hname,hname,MAXNPV,0,MAXNPV);

  //Get Areas//

  double geo[nEta][nEta];
  getGeometry(geo);

  //Set Branches//

  float weight = 1.;

  if (ISMC && REWEIGHT)
    T->SetBranchAddress("weight", &weight);

  ULong64_t event;
  int run, lumi, bx;

  if (!ISMC){
    //parsePileUpJSON2();
    //T->SetBranchAddress("event", &event);
    //T->SetBranchAddress("run", &run);
    //T->SetBranchAddress("lumi", &lumi);
    //T->SetBranchAddress("bx", &bx);
  }

  UChar_t f[numFlavors][nEta];
  float energy[nEta], eRMS[nEta], et[nEta];
  float mu, rho, rhoC0, rhoCC;
  int nPVall, nPV;
  float pv_ndof[nPVall], pv_z[nPVall], pv_rho[nPVall];
  float ht;
  int nJets=4;
  float jet_eta[nJets], jet_phi[nJets], jet_pt[nJets], jet_area[nJets];
  float jet_ch[nJets], jet_nh[nJets], jet_ne[nJets], jet_hfh[nJets], jet_hfe[nJets], jet_lep[nJets];

  T->SetBranchAddress("energy", energy);
  //T->SetBranchAddress("eRMS", eRMS);
  //T->SetBranchAddress("et", et);
  T->SetBranchAddress("fchm", f[chm]);
  T->SetBranchAddress("fchu", f[chu]);
  T->SetBranchAddress("fnh", f[nh]);
  T->SetBranchAddress("fne", f[ne]);
  T->SetBranchAddress("fhfh", f[hfh]);
  T->SetBranchAddress("fhfe", f[hfe]);
  T->SetBranchAddress("flep", f[lep]);
  T->SetBranchAddress("funtrk", f[untrk]);
  T->SetBranchAddress("mu", &mu);
  T->SetBranchAddress("rho", &rho);
  //T->SetBranchAddress("rhoC0", &rhoC0);
  //T->SetBranchAddress("rhoCC", &rhoCC);
  //T->SetBranchAddress("nPVall", &nPVall);
  T->SetBranchAddress("nPV", &nPV);  

  //Loop Over Entries//

  for (Long64_t n=0; n<nEntries; n++){
    if (n % 100000 == 0)
      cout << "Processing Event " << n+1 << endl;
    T->GetEntry(n);

    FillHist1D("nPU", mu, weight);
    FillHist1D("nPV", nPV, weight);
    FillHist1D("rho", rho, weight);
    FillProfile("p_nPV_nPU", mu, nPV, weight);
    FillProfile("p_rho_nPU", mu, rho, weight);
    FillProfile("p_rho_nPV", nPV, rho, weight);

    int intmu = mu + 0.5;

    for (int ieta=0; ieta<nEta; ieta++){
      double eta = 0.5*(etabins[ieta] + etabins[ieta+1]);

      for (int i_id=0; i_id<numFlavors; i_id++){
        double offpt = 0;

        for(int jeta = ieta-10; jeta <= ieta+10; jeta++){
          if( jeta<0 || jeta+1 >= nEta) continue;

          offpt += energy[jeta] * geo[ieta][jeta] * f[i_id][jeta] / 255.;
        }
        hname = Form("p_offset_eta_nPU%i_", intmu) + ids[i_id];
        FillProfile(hname, eta, offpt, weight);
        hname = Form("p_offset_eta_nPV%i_", nPV) + ids[i_id];
        FillProfile(hname, eta, offpt, weight);
      }
    }
  } //end event loop

  //Write Histos//

  TFile* outFile = new TFile(OUTFILENAME,"RECREATE");
  outFile->cd();

  for (int i_id=0; i_id<numFlavors; i_id++){
    outFile->mkdir( "offset_nPU/" + ids[i_id] );
    outFile->mkdir( "offset_nPV/" + ids[i_id] );
  }

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++)
    hid->second->Write();

  for (map<TString, TH2*>::iterator hid = m_Histos2D.begin(); hid != m_Histos2D.end(); hid++)
    hid->second->Write();

  for (map<TString, TProfile*>::iterator hid = m_Profiles.begin(); hid != m_Profiles.end(); hid++){
    outFile->cd();
    hname = hid->first;

    if ( hname.Contains("p_offset_eta_nPU") ){
      TString id = hname( hname.Last('_')+1, hname.Length() );
      outFile->cd(OUTFILENAME + ":/offset_nPU/" + id);
    }
    else if ( hname.Contains("p_offset_eta_nPV") ){
      TString id = hname( hname.Last('_')+1, hname.Length() );
      outFile->cd(OUTFILENAME + ":/offset_nPV/" + id);
    }
    hid->second->Write();
  }

  for (map<TString, TProfile2D*>::iterator hid = m_Profiles2D.begin(); hid != m_Profiles2D.end(); hid++)
    hid->second->Write();

  outFile->Write();
  delete outFile;
  outFile = 0;
}

//Functions//

void FillHist1D(const TString& histName, const Double_t& value, double weight) 
{
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, double weight) 
{
  map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    cout << "%FillHist2D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void FillProfile(const TString& histName, const Double_t& value1, const Double_t& value2, double weight) 
{
  map<TString, TProfile*>::iterator hid=m_Profiles.find(histName);
  if (hid==m_Profiles.end())
    cout << "%FillProfile -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void FillProfile2D(const TString& histName, const Double_t& value1, const Double_t& value2, const Double_t& value3, double weight) 
{
  map<TString, TProfile2D*>::iterator hid=m_Profiles2D.find(histName);
  if (hid==m_Profiles2D.end())
    cout << "%FillProfile2D -- Could not find profile with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, value3, weight);
}

void getGeometry(double (&geo)[nEta][nEta]){

  for (int ieta=0; ieta<nEta; ieta++){
    double eta = 0.5*(etabins[ieta] + etabins[ieta+1]);

    for(int jeta = ieta-10; jeta <= ieta+10; jeta++){   

      if( jeta<0 || jeta+1 >= nEta) continue;
    
      double etaL = etabins[jeta];                    // left  edge of the eta strip
      double etaR = etabins[jeta+1];                  // right edge of the eta strip
      double etaC = 0.5*(etaL+etaR);                  // center of the eta strip
      double A = areaS(rCone, etaL-eta, etaR-eta);    // area of the eta strip inside the cone
      //if (A <= 0.00000001) continue;                //  (A <= 0.)  gives precision issues

      // the next lines make sure that we do vector addition in phi direction;
      // We would be doing scalar addition with coef = 1.

      double dphi = dist(rCone, etaL-eta, etaR-eta);
      double coef = dphi > 0.000001 ? TMath::Sin(dphi) / dphi : 1.;

      geo[ieta][jeta] = A * coef / (2*TMath::Pi()*(etaR-etaL)) / cosh(etaC);
    }
  }
}

double areaS(double R, double x1, double x2){
//
// Area inside a shape delinated by a circle of radius R
// centered at (0,0) and two parallel lines at x=x1 and x=x2
//

   if( R<=0. || x1==x2 ) return 0.;
   if(x1*x2>0 && fabs(x1) > R && fabs(x2) > R) return 0. ;

   double d1 = fabs(x1) ;
   double d2 = fabs(x2) ;
   if (d1>d2){
     double d = d1;
     d1 = d2;
     d2 = d ;
   }

   // area of segment at distance d1 from the center
   double theta1 = 2.* TMath::ACos(d1/R) ;
   double A1 = (0.5*R*R*(theta1-TMath::Sin(theta1))) ;

   // area of segment at distance d2 from the center
   double A2 = 0. ;
   if(d2<=R){
     double theta2 = 2.* TMath::ACos(d2/R) ;
     A2 = (0.5*R*R*(theta2-TMath::Sin(theta2))) ;
   }

   if(x1*x2>=0){ // both lines on the same side from the center
     return A1 - A2 ;
   } else { // the lines on the opposite side from the center
    return TMath::Pi()*R*R - A1 - A2 ;
   }
}

double dist(double R, double x1, double x2){
//
// Take a circle of radius R centered at (0,0).
// This function calculates distance between a 
// midpoint of x=x1 and x=x2 and the point on circle rim 
// along vertical line.
//

   if( R<=0.) return 0.;
   if(x1*x2>0 && fabs(x1) >= R && fabs(x2) >= R) return 0. ; // both lines outside the circle 
                                                             // and on the same side of origin. 
  if(fabs(x1)>R) x1 = TMath::Sign(R, x1);
  if(fabs(x2)>R) x2 = TMath::Sign(R, x2); 

   double x = 0.5*( x1 + x2) ;

   return TMath::Sqrt(R*R-x*x);
}
