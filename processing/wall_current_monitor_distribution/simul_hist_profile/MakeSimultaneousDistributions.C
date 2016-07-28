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

  TFile* out_root = new TFile(output_file.c_str(),"RECREATE");

  // Load files to list
  std::ifstream in_file(input_list.c_str());
  if(!in_file){
    std::cout << "couldn't open " << input_list << std::endl;
    return 1;
  }
  std::string line = "";

  // Loop over profiles for same run to adjust them simultaneously
  while(getline(in_file,line)){
    std::string blue_profile_dat,yell_profile_dat,blue_profile_name,yell_profile_name;
    std::stringstream ss;
    ss.str(line.c_str());
    ss >> blue_profile_dat >> blue_profile_name >> yell_profile_dat >> yell_profile_name;
//    std::cout <<
//      blue_profile_dat << std::endl 
//      << blue_profile_name << std::endl 
//      << yell_profile_dat << std::endl
//      << yell_profile_name << std::endl << std::endl;
    
    // Blue Profile
    std::string bl_distribution = blue_profile_dat;
    std::string bl_profile_name = blue_profile_name;
    std::string bl_profile_name_tail = bl_profile_name + "_tail";

    //std::cout << bl_distribution << ", " << bl_profile_name << std::endl;
    std::vector<double>x_bl;
    std::vector<double>y_bl;
    LoadWCMProfile(bl_distribution,x_bl,y_bl);
    float bin_w_b = x_bl[1] - x_bl[0];
    
    std::string ye_distribution = yell_profile_dat;
    std::string ye_profile_name = yell_profile_name;
    std::string ye_profile_name_tail = ye_profile_name + "_tail";

    //std::cout << ye_distribution << ", " << ye_profile_name << std::endl;
    std::vector<double>x_ye;
    std::vector<double>y_ye;
    LoadWCMProfile(ye_distribution,x_ye,y_ye);
    float bin_w_y = x_ye[1] - x_ye[0];

    auto ye_y_max = std::max_element(y_ye.begin(),y_ye.end());
    auto bl_y_max = std::max_element(y_bl.begin(),y_bl.end());
    unsigned int max_pos_blue = std::distance(y_bl.begin(),bl_y_max);
    unsigned int max_pos_yell = std::distance(y_ye.begin(),ye_y_max);

    std::cout << bin_w_y << ", " << bin_w_b << ", " << fabs(bin_w_b - bin_w_y) << std::endl;
    std::cout << y_bl.size() << ", " << y_ye.size() << std::endl;
    std::cout << x_bl[max_pos_blue] << ", " << x_ye[max_pos_yell] << ", " << fabs(x_bl[max_pos_blue] - x_ye[max_pos_yell]) <<  std::endl;

    // To do: overlap at the average bin between the max value of the two
    // profiles, build out profiles to either side of this. Fill to histograms
    // for visualization, but actually save to data files, intripolate as
    // normal.

  }
  out_root->Write();
  out_root->Close();
  delete out_root;
  return 0;
}
