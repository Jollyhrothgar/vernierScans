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
    
    // Set up the names
    std::string bl_distribution = blue_profile_dat;
    std::string bl_profile_name = blue_profile_name;
    std::string ye_distribution = yell_profile_dat;
    std::string ye_profile_name = yell_profile_name;

    // Store the data
    std::vector<double>x_bl;
    std::vector<double>y_bl;
    std::vector<double>x_ye;
    std::vector<double>y_ye;
    LoadWCMProfile(bl_distribution,x_bl,y_bl);
    LoadWCMProfile(ye_distribution,x_ye,y_ye);

    // Obtain the bin width of the data
    float bin_w_b = x_bl[1] - x_bl[0];
    float bin_w_y = x_ye[1] - x_ye[0];

    // Get the nominal zero point by taking the average 
    // between the max of the two distributions. The overlap point is set
    // to be the midway bin between.
    auto ye_y_max = std::max_element(y_ye.begin(),y_ye.end());
    auto bl_y_max = std::max_element(y_bl.begin(),y_bl.end());
    unsigned int max_pos_blue = std::distance(y_bl.begin(),bl_y_max);
    unsigned int max_pos_yell = std::distance(y_ye.begin(),ye_y_max);

    unsigned int center = abs(max_pos_blue + max_pos_yell)/2;
    float step = bin_w_y;

    // Shift the profiles to lie at z = 0, make them have a range of +/- 6000
    std::deque<float> blue_profile_x;
    std::deque<float> blue_profile_y;
    std::deque<float> yell_profile_x;
    std::deque<float> yell_profile_y;

    blue_profile_x.push_front(0.0);
    yell_profile_x.push_front(0.0);
    blue_profile_y.push_front(y_bl[center]);
    yell_profile_y.push_front(y_ye[center]);

    // Fill in the left
    auto blue_itr = y_bl.begin()+center-1;
    auto yell_itr = y_ye.begin()+center-1;
    float x = 0. - step;
    while(blue_itr != y_bl.begin() && yell_itr != y_ye.begin()){
      blue_profile_x.push_front(x);
      yell_profile_x.push_front(x);
      blue_profile_y.push_front(fabs(*blue_itr));
      yell_profile_y.push_front(fabs(*yell_itr));

      x-=step;
      if( fabs(x) > 450) {
        break;
      }
      --blue_itr;
      --yell_itr;
    }
    while(x > -10000.0){
      blue_profile_x.push_front(x);
      yell_profile_x.push_front(x);
      blue_profile_y.push_front(0.0);
      yell_profile_y.push_front(0.0);
      x-=step;
    }

    // Fill in the right
    x = 0.+step;
    blue_itr = y_bl.begin()+center+1;
    yell_itr = y_ye.begin()+center+1;
    while(blue_itr != y_bl.end() && yell_itr != y_ye.end()){
      blue_profile_x.push_back(x);
      yell_profile_x.push_back(x);
      blue_profile_y.push_back(fabs(*blue_itr));
      yell_profile_y.push_back(fabs(*yell_itr));
      x+=step;
      
      if( fabs(x) > 450) {
        break;
      }
      ++blue_itr;
      ++yell_itr;
    }
    while(x < 10000.0){
      blue_profile_x.push_back(x);
      yell_profile_x.push_back(x);
      blue_profile_y.push_back(0.0);
      yell_profile_y.push_back(0.0);
      x+=step;
    }
    
    TH1F* blue_profile = new TH1F(bl_profile_name.c_str(),
        bl_profile_name.c_str(),
        blue_profile_x.size(),
        *std::min_element(blue_profile_x.begin(),blue_profile_x.end())-step,
        *std::max_element(blue_profile_x.begin(),blue_profile_x.end())+step
        );
    TH1F* yell_profile = new TH1F(ye_profile_name.c_str(),
        ye_profile_name.c_str(),
        yell_profile_x.size(),
        *std::min_element(yell_profile_x.begin(),yell_profile_x.end())-step,
        *std::max_element(yell_profile_x.begin(),yell_profile_x.end())+step
        );
    for(unsigned int bin_i = 0; bin_i < blue_profile_y.size(); bin_i++){
      blue_profile->SetBinContent(bin_i+1,blue_profile_y[bin_i]);
      yell_profile->SetBinContent(bin_i+1,yell_profile_y[bin_i]);
    }
    
    blue_profile->Smooth(300);
    blue_profile->Scale(1.0/blue_profile->Integral("width"));
    yell_profile->Smooth(300);
    yell_profile->Scale(1.0/yell_profile->Integral("width"));

    blue_profile->Write();
    yell_profile->Write();
    delete blue_profile;
    delete yell_profile;
  }
  out_root->Write();
  out_root->Close();
  delete out_root;
  return 0;
}
