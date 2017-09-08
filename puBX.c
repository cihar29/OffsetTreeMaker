void puBX() {

  ifstream file("pu_lumi_bx.txt");
  ofstream outfile("short.txt");
  string line;

  while (getline(file, line)) {
    if (line.at(0) == '#') continue;

    string str, run, ls, bxs = "";
    int delim_pos;

    for (int string_num=0; (delim_pos = line.find(",")) != -1; string_num++) {

      str = line.substr(0, delim_pos);
      line.erase(0, delim_pos + 1);

      if (string_num == 0)      run = str.substr(0, str.find(":"));
      else if (string_num == 1) ls = str.substr(0, str.find(":"));
      else if (string_num == 8) {
        line.erase(0, 1);  //get rid of '['
        line.erase(line.size()-1, 1);
        line += " ";       //replace ']' with space

        while ( (delim_pos = line.find(' ')) != -1) {
          string id = line.substr(0, delim_pos);
          line.erase(0, delim_pos + 1);

          delim_pos = line.find(' ');    //delivered
          line.erase(0, delim_pos + 1);

          delim_pos = line.find(' ');    //recorded
          string lumi = line.substr(0, delim_pos);

          if (stof(lumi) > 1.) bxs += id + " " + lumi + " ";
          line.erase(0, delim_pos + 1);
        }
      }
    }
    outfile << run << " " << ls << " " << bxs << endl;
  }
  file.close();
  outfile.close();
}
