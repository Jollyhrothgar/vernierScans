#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"

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
    std::vector<double> distance_cm_;
    std::vector<double> wcm_intensity_;
    std::vector<double> noise;

    std::string data = "";
    while(getline(in_data,data)) {
      std::stringstream wcm_ss;
      wcm_ss.str(data);
      double nanoseconds_;
      double intensity;
      double distance;

      wcm_ss >> nanoseconds_ >> intensity;

      // Transform
      // nanoseconds_ = nanoseconds_*-1.0+105.14;
      wcm_intensity_.push_back(intensity);
      time_ns_.push_back(nanoseconds_);
      distance_cm_.push_back(nanoseconds_*to_seconds*c_vel);

      if(nanoseconds_ > time_cut) {
        noise.push_back(intensity);
      }
    }
    // Goal: create real model of WCM profile which includes the real
    // fluctuation of the beam gas which could impact the data. Though this
    // affect is small, it is relatively easy to implement. 
    
    // GENERATE CORRECTED WCM PROFILE, ALL POSITIVE
    std::stringstream name;
    std::stringstream title;
    std::vector<double> density;
    std::vector<double> _time = time_ns_;
    double offset = *std::min_element(wcm_intensity_.begin(), wcm_intensity_.end());
    
    for(unsigned int entry_i = 0; entry_i < wcm_intensity_.size(); entry_i++) {
      density.push_back(wcm_intensity_[entry_i] + fabs(offset));
    }
    name.str("");
    title.str("");
    TGraph* g_wcm_norm = new TGraph(time_ns_.size(),&time_ns_[0],&density[0]);
    name << dist_name << "_time_shifted";
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
      *std::min_element(noise.begin(),noise.end())+fabs(offset),
      *std::max_element(noise.begin(),noise.end())+fabs(offset)
      );

    for(auto noise_i = noise.begin(); noise_i != noise.end(); ++noise_i) {
      h_fluctuation->Fill(*noise_i + fabs(offset));
    }
    h_fluctuation->Smooth(3);
    h.push_back(h_fluctuation);

    // GENERATE CENTERED WCM PROFILE, WITH NOISE SAMPLING TO FILL EMPTY BUNCH

    // 1. Split profile from the rest
    std::vector<double> z_prof;
    std::vector<double> left_noise;
    std::vector<double> right_noise;
    for(unsigned int entry_i = 0; entry_i < time_ns_.size(); entry_i++ ) {
      if(time_ns_[entry_i] < time_cut) {
        z_prof.push_back(density[entry_i]);
	// Get bounding distributions of random noise, sampled from the real
	// noise of the WCM profile (hypothesis: this isn't noise, its beam gas,
	// so instead, lets just shift the WCM profile, keeping as much data as
	// possible. The beam is a ring anyway, so periodic boundary conditions
	// make sense. This is kept in case cross-checking is ever desired.
        // left_noise.push_back(h_fluctuation->GetRandom());
        // right_noise.push_back(h_fluctuation->GetRandom());
      } else if((time_ns_[entry_i] >= time_cut)&&(time_ns_[entry_i] < time_cut*2.0)) {
        right_noise.push_back(density[entry_i]);
      } else {
        left_noise.push_back(density[entry_i]);
      }
    }

    // Now, create the profile
    std::vector<double> density_shifted;
    density_shifted.insert(density_shifted.end(), left_noise.begin(), left_noise.end());
    density_shifted.insert(density_shifted.end(), z_prof.begin(), z_prof.end());
    density_shifted.insert(density_shifted.end(), right_noise.begin(), right_noise.end());

    auto density_itr = density_shifted.begin();
    auto time_itr = time_ns_.begin();
    name.str("");
    title.str("");
    TGraph* z_bunch_prof = new TGraph();
    name << dist_name << "_bunch_profile";
    title << dist_name << ";Z-position (PHENIX IR);density";
    z_bunch_prof->SetName(name.str().c_str());
    z_bunch_prof->SetTitle(title.str().c_str());
    while( (density_itr != density_shifted.end() ) && ( time_itr != time_ns_.end()) ) {
      double time_shifted = *time_itr;
      time_shifted = time_shifted - ((time_cut*3.0)/2.0); // Shift so bunch is centered at z = 0
      z_bunch_prof->SetPoint(z_bunch_prof->GetN(), time_shifted*to_seconds*c_vel, *density_itr);
      ++density_itr;
      ++time_itr;
    }
    g.push_back(z_bunch_prof);

    // Now, create the normalized profile
    //
    // Create a textfile for alternative usage
    //
    name.str("");
    title.str("");
    TGraph* z_prof_norm = new TGraph();
    name << dist_name << "_density";
    title << dist_name << ";Z-position (PHENIX IR) cm;P.D.F. Density";
    z_prof_norm->SetName(name.str().c_str());
    z_prof_norm->SetTitle(title.str().c_str());
    double norm = z_bunch_prof->Integral();
    density_itr = density_shifted.begin();
    time_itr = time_ns_.begin();

    name.str("");
    name << dist_name << "_density.txt";

    std::ofstream wcm_out(name.str().c_str()); 
    while( (density_itr != density_shifted.end() ) && ( time_itr != time_ns_.end()) ) {
      double time_shifted = *time_itr;
      time_shifted = time_shifted - ((time_cut*3.0)/2.0); // Shift so bunch is centered at z = 0
      z_prof_norm->SetPoint(z_prof_norm->GetN(), (time_shifted*to_seconds*c_vel), *density_itr/norm);
      wcm_out << (time_shifted*to_seconds*c_vel) << " " << *density_itr/norm << std::endl;
      ++density_itr;
      ++time_itr;
    }
    wcm_out.close();
    g.push_back(z_prof_norm);

    name.str("");
    title.str("");
    TGraph* g_wcm = new TGraph(time_ns_.size(),&time_ns_[0],&wcm_intensity_[0]);
    name << dist_name << "_time_raw";
    title << name.str() << ";nanoseconds;intensity";
    g_wcm->SetName(name.str().c_str());
    g_wcm->SetTitle(title.str().c_str());
    g.push_back(g_wcm);

    g_wcm = new TGraph(time_ns_.size(), &distance_cm_[0], &wcm_intensity_[0]);
    name.str("");
    title.str("");
    name << dist_name << "_space_raw";
    title << name.str() << ";centimeters;intensity";
    g_wcm->SetName(name.str().c_str());
    g_wcm->SetTitle(title.str().c_str());
    g.push_back(g_wcm);
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
