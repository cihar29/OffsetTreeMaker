#ifndef __parsePileUpJSON2_C__
#define __parsePileUpJSON2_C__

#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace std;

map<int, map<int, map<int, float> > > m_PU;  //m_PU[run][ls][bx]
const float MINBIAS_XS = 0.692;

float getAvgPU(int run, int ls, int bx) {

  return m_PU[run][ls][bx];
}

int parsePileUpJSON2(string filename="short_pileup_9_7_2016.txt") {

  //### Using Brilcalc ###//
  cout << "Opening " << filename << "...";

  string line;
  ifstream file(filename);

  if (file.is_open()) {
    cout << "ok" << endl;

    //loop over lines in file
    while ( getline(file,line) ) {

      int delim_pos = line.find(' ');
      if (delim_pos == -1) continue;

      string run_str = line.substr(0, delim_pos);
      line.erase(0, delim_pos + 1);

      delim_pos = line.find(' ');
      string ls_str = line.substr(0, delim_pos);
      line.erase(0, delim_pos + 1);

      map<int, float> bx_PU;

      while ( (delim_pos = line.find(' ')) != -1) {
        string bx = line.substr(0, delim_pos);
        line.erase(0, delim_pos + 1);

        delim_pos = line.find(' ');
        string lumi = line.substr(0, delim_pos);
        line.erase(0, delim_pos + 1);

        bx_PU[stoi(bx)] = stof(lumi) * MINBIAS_XS;
      }
      m_PU[ stoi(run_str) ][ stoi(ls_str) ] = bx_PU;
    }
    file.close();
  }
  else
    cout << "Unable to open file" << endl;

  return 0;
}

#endif //__parsePileUpJSON2_C__
