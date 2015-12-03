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

int main( int argc, char** argv ) {

  // Give a list of WCM profiles to plot
  std::string input_list = argv[1];
  std::string output_file = argv[2];

  if(argc != 3) {
    std::cout << "usage is: " << argv[0] << " <list_file> <output_file>" << std::endl;
    return 1;
  }

  const double time_cut = 35.0467;
  const double c_vel = 2.998e+10; // speed of light in cm/s
  double to_seconds = 1.0e-9;     // entry in data file is in nanoseconds
  std::ifstream in_list(input_list.c_str());

  std::vector<TGraph*> g;// all the graphs to write out
  std::vector<TH1F*> h; // all the histograms to write out

  // Slice and translate, filling with random noise

  if( !in_list ) return 1;
  std::string line = "";
  while(getline(in_list,line)) {
    std::stringstream ss;
    ss.str(line);
    std::string wcm_data;
    std::string dist_name;
    ss >> wcm_data >> dist_name;

    std::ifstream in_data(wcm_data.c_str());
    if(!in_data) return 1;

    // reset each time (because it leaves scope)
    std::vector<double> time_ns_;
    std::vector<double> wcm_intensity_;
    std::vector<double> noise;

    std::string data = "";
    while(getline(in_data,data)) {
      std::stringstream wcm_ss;
      wcm_ss.str(data);
      double nanoseconds_;
      double intensity;

      wcm_ss >> nanoseconds_ >> intensity;

      wcm_intensity_.push_back(intensity);
      time_ns_.push_back(nanoseconds_);

      if(nanoseconds_ > time_cut) {
        noise.push_back(intensity);
      }
    }
    
    // GENERATE CORRECTED WCM PROFILE, ALL POSITIVE
    std::stringstream name;
    std::stringstream title;
    std::vector<double> density;
    std::vector<double> _time = time_ns_;
    double offset = *std::min_element(wcm_intensity_.begin(), wcm_intensity_.end());
    
    for(unsigned int entry_i = 0; entry_i < wcm_intensity_.size(); entry_i++) {
      density.push_back(wcm_intensity_[entry_i] + fabs(offset));
    }

    for(unsigned int noise_i = 0; noise_i < noise.size(); noise_i++) {
      noise[noise_i] = noise[noise_i] + fabs(offset);
    }

    // Now we have our shifted intensity (called density) and our noise distribution (called noise)

    name.str("");
    title.str("");
    TGraph* g_wcm_norm = new TGraph(time_ns_.size(),&time_ns_[0],&density[0]);
    name << dist_name << "_time_pos_values";
    title << name.str() << ";nanoseconds;density";
    g_wcm_norm->SetName(name.str().c_str());
    g_wcm_norm->SetTitle(title.str().c_str());
    g.push_back(g_wcm_norm);

    // GENERATE NOISE DISTRIBUTION FOR WCM FILL
    name.str("");
    title.str("");
    name << dist_name << "_fluctuation";
    title << name.str() << ";density;counts";
    TH1F* h_fluctuation = new TH1F(
      name.str().c_str(),
      title.str().c_str(),
      100,
      *std::min_element(noise.begin(),noise.end()),
      *std::max_element(noise.begin(),noise.end())
      );

    for(auto noise_i = noise.begin(); noise_i != noise.end(); ++noise_i) {
      h_fluctuation->Fill(*noise_i);
    }
    h_fluctuation->Smooth(3);
    h.push_back(h_fluctuation);
    
    // GENERATE CENTERED WCM PROFILE, WITH NOISE SAMPLING TO FILL EMPTY BUNCH
    // 1. Split profile from the rest
    // 1.a. Find maximum of profile
    // 1.b. Slice around maximum
    // 1.c. expontential decay to 0 after main bunch edge
    auto max_entry = std::max_element(density.begin(),density.end());

    // time cut == 1/3 of the domain.
    double half_time_cut = time_cut/2.0;
    auto first_time = time_ns_.begin();
    auto second_time = time_ns_.begin();
    ++second_time;
    double time_interval = *second_time - *first_time;
    std::deque<double> density_centered;
    std::deque<double> time_centered;
    
    density_centered.push_front(*max_entry); // now the maximum value of the WCM is in the middle
    time_centered.push_front(0.0); // we set the time in the middle to be "0"
    double left_time = 0.;
    double right_time = 0.;
    left_time -= time_interval;
    right_time += time_interval;

    auto left_itr = max_entry;
    auto right_itr = max_entry;
    --left_itr;
    ++right_itr;

    name.str("");
    title.str("");
    name << dist_name << "_left_side";
    title << dist_name << " left side of bunch";
    TGraph* left_side = new TGraph();
    left_side->SetName(name.str().c_str());
    left_side->SetTitle(name.str().c_str());
    g.push_back(left_side);

    name.str("");
    title.str("");
    name << dist_name << "_right_side";
    title << dist_name << " right side of bunch";
    TGraph* right_side = new TGraph();
    right_side->SetName(name.str().c_str());
    right_side->SetTitle(name.str().c_str());
    g.push_back(right_side);

    while( left_time > (-1.*half_time_cut) && right_time < half_time_cut ) {
      density_centered.push_front(*left_itr);
      density_centered.push_back(*right_itr);
      time_centered.push_front(left_time);
      time_centered.push_back(right_time);
      --left_itr;
      ++right_itr;
      left_time -= time_interval;
      right_time += time_interval;

      // if true, we are in tails of distribution
      if(*left_itr  < *max_entry/5.0) {
        left_side->SetPoint(left_side->GetN(),left_time,*left_itr);
      } 
      if(*right_itr < *max_entry/5.0) {
        right_side->SetPoint(right_side->GetN(),right_time,*right_itr);
      } 
      if(left_itr == density.begin()) break; // left side will hit the edge before the right side
    }
    while( left_time > (-1.0*(time_cut*3./2.)) && right_time < (time_cut*3./2.) ) {
      density_centered.push_front(0.0);
      density_centered.push_back (0.0);
      time_centered.push_front(left_time);
      time_centered.push_back(right_time);
      left_time -= time_interval;
      right_time += time_interval;
    }

    name.str("");
    title.str("");
    TGraph* centered_density_ = new TGraph();
    name << dist_name << "_max_centered";
    title << dist_name << " centered on maximum;position (cm);WCM Intensity";
    centered_density_->SetName(name.str().c_str());
    centered_density_->SetTitle(name.str().c_str());
    auto d_i = density_centered.begin();
    auto t_i = time_centered.begin();
    while( d_i != density_centered.end() && t_i != time_centered.end() ) {
      centered_density_->SetPoint(centered_density_->GetN(), ToCentimeters(*t_i), fabs(*d_i));
      ++d_i;
      ++t_i;
    }
    g.push_back(centered_density_);

    double norm = centered_density_->Integral();
    std::cout << "Normalization: " << norm << std::endl;

    // Now, create the normalized profile
    // Create a textfile for alternative usage
    name.str("");
    title.str("");
    TGraph* z_prof_norm = new TGraph();
    name << dist_name << "_density";
    title << dist_name << ";Z-position (PHENIX IR) cm;P.D.F. Density";
    z_prof_norm->SetName(name.str().c_str());
    z_prof_norm->SetTitle(title.str().c_str());
    
    auto density_c_itr = density_centered.begin();
    auto time_c_itr = time_centered.begin();

    name.str("");
    name << "./profile_densities/" << dist_name << "_density.txt";

    std::ofstream wcm_out(name.str().c_str()); 
    while( (density_c_itr != density_centered.end() ) && ( time_c_itr != time_centered.end()) ) {
      z_prof_norm->SetPoint(z_prof_norm->GetN(), ToCentimeters(*time_c_itr), fabs(*density_c_itr/norm));
      wcm_out << ToCentimeters(*time_c_itr) << " " << fabs(*density_c_itr/norm) << std::endl;
      ++density_c_itr;
      ++time_c_itr;
    }
    wcm_out.close();
    g.push_back(z_prof_norm);
  } // end loop over files

  TFile* f = new TFile(output_file.c_str(), "RECREATE");
  
  for(auto g_itr = g.begin(); g_itr != g.end(); ++g_itr) {
    (*g_itr)->Write();
    delete (*g_itr);
  }
  for(auto h_itr = h.begin(); h_itr != h.end(); ++h_itr) {
    (*h_itr)->Write();
    delete(*h_itr);
  }
  f->Write();
  delete f;
  return 0;
}
