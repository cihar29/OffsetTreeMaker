#ifndef __parsePileUpJSON2_C__
#define __parsePileUpJSON2_C__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"

using namespace std;

map<int, map<int, float> > m_PU;
map<int, map<int, map<int, float> > > m_PUbx;
const float MINBIAS_XS = 69200;

float getAvgPU(int run, int ls) {
  return m_PU[run][ls];
}

//float getAvgPU(int run, int ls, int bx) {
//  return m_PUbx[run][ls][bx];
//}
/*
int parsePileUpJSON2(string filename="lumi-per-bx.root") {

  //### Lumi per BX ###//

  TFile* file = TFile::Open( filename.data() );
  int nBuckets = 3564;

  TIter nextkey(file->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey*)nextkey()) ) {
    TString keyname = key->GetName();
    TTree* tree = (TTree*) file->Get(keyname);
    int run = stoi( keyname.Data() );

    UInt_t ls;
    float lumi[nBuckets];

    tree->SetBranchAddress("LumiSection", &ls);
    tree->SetBranchAddress("Luminosity",  lumi);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t n=0; n<nEntries; n++) {
      tree->GetEntry(n);

      for (int i=0; i<nBuckets; i++) m_PUbx[run][ls][i] = lumi[i] * MINBIAS_XS / pow(2,18);
    }
  }

*/
int parsePileUpJSON2(string filename="pileup.txt") {
  //### Using Brilcalc ###//
/*  cout << "Opening " << filename << "...";

  string line;
  ifstream file(filename);

  if (file.is_open()){
    cout << "ok" << endl;

    //loop over lines in file
    while ( getline(file,line) ){

      string run_str, ls_str;
      string str;
      int delim_pos;
      float PU = -1;

      if ( line.at(0) != '#' ){

        //loop over strings in line
        for (int string_num=0; (delim_pos = line.find(",")) != -1; string_num++){

          str = line.substr(0, delim_pos);
          line.erase(0, delim_pos + 1);

          if (string_num == 0)  //first string holds run number
            run_str = str.substr(0, str.find(":"));

          else if (string_num == 1) //second string has ls
            ls_str = str.substr(0, str.find(":"));

          else if (string_num == 7) //eighth string has pu
            PU = stof( str );
        }
        int run = stoi( run_str );
        int ls = stoi( ls_str );

        m_PU[run][ls] = PU;
      }
    }
    file.close();
  }
  else
    cout << "Unable to open file" << endl;

*/
  //### Using Pileup JSON file ###//
  cout << "Minimum Bias Cross Section: " << MINBIAS_XS << endl;
  cout << "Opening " << filename << "...";

  string line;
  ifstream file(filename);

  if (file.is_open()){
    cout << "ok" << endl;

    //loop over lines in file
    while ( getline(file,line) ){

      string str;
      int delim_pos, run, ls;
      float PU = -1;

      if ( line.at(0) != '#' ){

        //loop over strings in line
        for (int string_num=0; (delim_pos = line.find(" ")) != -1; string_num++){

          str = line.substr(0, delim_pos);
          line.erase(0, delim_pos + 1);

          if (string_num == 0)  //first string holds run number
            run = stoi( str );

          else if (string_num == 1) //second string has ls
            ls = stoi( str );
        }
        PU = stod( line ) * MINBIAS_XS;
        m_PU[run][ls] = PU;
      }
    }
    file.close();
  }
  else
    cout << "Unable to open file" << endl;

  return 0;
}

#endif //__parsePileUpJSON2_C__

