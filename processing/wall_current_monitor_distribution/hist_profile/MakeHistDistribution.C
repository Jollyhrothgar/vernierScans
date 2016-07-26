#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <deque>
#include <iterator>

#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"

double ToCentimeters( double nanoseconds ) {
  return nanoseconds*1.0e-9*2.998e+10;
}

int LoadFileList(const std::string& list, std::vector<std::pair<std::string,std::string> > &v ){
  std::ifstream in_list(list.c_str());
  std::string line;
  if (!in_list) {
    std::cout << "couldn't find file, " << list << std::endl;
    return 0;
  }
  while(getline(in_list,line)){
    std::string filename;
    std::string outname;
    std::stringstream ss;
    ss.str(line);
    ss >> filename >> outname;
    std::cout << filename << ", " << outname << std::endl;
    v.push_back(std::make_pair(filename,outname));
  }
  return 1;
}

int LoadWCMProfile(const std::string& name, std::vector<double>&x,std::vector<double>&y){
  std::ifstream in_file(name.c_str());
  if (!in_file) {
    std::cout << "couldn't find file, " << name << std::endl;
    return 0;
  }
  std::string line;
  while(getline(in_file,line)){
    std::stringstream ss;
    ss.str(line);
    double xval,yval;
    ss >> xval >> yval;
    x.push_back(ToCentimeters(xval));
    y.push_back(yval);
  }
  return 0;
}

int main( int argc, char** argv ) {
  // Give a list of WCM profiles to plot
  if(argc != 3) {
    std::cout << "usage is: " << argv[0] << " <list_file> <output_file>" << std::endl;
    return 1;
  }
  std::string input_list = argv[1];
  std::string output_file = argv[2];

  // Load files to list
  std::vector<std::pair<std::string,std::string> > files;
  LoadFileList(input_list,files);

  TFile* out_root = new TFile(output_file.c_str(),"RECREATE");

  for(auto file = files.begin(); file != files.end(); ++file){
    std::string in_distribution = file->first;
    std::string profile_name = file->second;
    std::string profile_name_tail = profile_name + "_tail";

    std::cout << in_distribution << ", " << profile_name << std::endl;
    std::vector<double>x;
    std::vector<double>y;
    LoadWCMProfile(in_distribution,x,y);
    float bin_w = x[1] - x[0];
    TH1F* profile = new TH1F(
        profile_name.c_str(),
        profile_name.c_str(),
        x.size(),
        *std::min_element(x.begin(),x.end())-0.5*bin_w,
        *std::max_element(x.begin(),x.end())+0.5*bin_w
        );
    TH1F* profile_tail = new TH1F(
        profile_name_tail.c_str(),
        profile_name_tail.c_str(),
        100,-0.1,0.1
        );
    for(unsigned int element = 0; element < x.size(); element++){
      if (x[element] > 1500.){
        profile_tail->Fill(y[element]);
      }
    }
    for(unsigned int element = 0; element < x.size(); element++){
      double offset = fabs(profile_tail->GetMean());
      double profile_value = y[element]+offset;
      profile_value = fabs(profile_value);
      profile->SetBinContent(element+1,profile_value);
    }
    out_root->cd();
    profile->Rebin(5);
    profile->Smooth(100);
    profile->Scale(1.0/profile->Integral("width"));
    profile->Write();
    profile_tail->Write();
    delete profile;
    delete profile_tail;
  }
  out_root->Write();
  out_root->Close();
  delete out_root;
  return 0;
}
