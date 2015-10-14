#include "HourglassData.h"
#include "MergedPRDFDST.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>

#include "TFile.h"
#include "TCanvas.h"

HourglassData::HourglassData() {
  std::cout << "Instantiating instance of HourglassData at " << this << std::endl;
}

HourglassData::~HourglassData() {
  std::cout << "Destroying instance of HourglassData at: " << this << std::endl;
  for(auto i = plot_registry_.begin(); i != plot_registry_.end(); ++i) {
    auto deleteme = (*i);
    if(deleteme) {
      delete deleteme;
    }
  }
}

int HourglassData::Init(
  const std::string& run_number,
  const std::string& scalers_file,
  const std::string& epoch_steps_file,
  const std::string& bpm_data_file_name, 
  const std::string& bpm_steps_file_name
){
  std::cout << "Processing Hourglass effect for run: " << run_number << std::endl;
  run_number_ = run_number;
  std::cout << " Reading infomration from " << scalers_file << " and " << std::endl
      << "    " << epoch_steps_file << std::endl;

  epoch_steps_file_name_ = epoch_steps_file;
  scalers_file_name_     = scalers_file;
  bpm_.Init(
    run_number_,
    bpm_data_file_name,
    bpm_steps_file_name
     );
  bpm_.Run();
  // Get micron position associated with step
  double h_1 = 0.;
  double h_2 = 0.;
  double v_1 = 0.;
  double v_2 = 0.;
  for(unsigned int i = 0; i < bpm_.GetNumberOfSteps(); i++) {
    h_steps_.push_back(bpm_.GetHStep(i));
    v_steps_.push_back(bpm_.GetVStep(i));
    if(i < bpm_.GetNumberOfSteps() / 2) { // first half
      h_1 += fabs(h_steps_.back());
      v_1 += fabs(v_steps_.back());
    } else { // second half
      h_2 += fabs(h_steps_.back());
      v_2 += fabs(v_steps_.back());
    }
  }
  bool horizontal_first = false;
  if( h_1 > h_2 ) { // horiontal scan first
    horizontal_first = true; 
  } else { // vertical scan first
    horizontal_first = false;
  }
  for (unsigned int i = 0; i < bpm_.GetNumberOfSteps(); i++) {
    double displacement = 0;
    displacement = pow(pow(h_steps_[i],2.0)+pow(v_steps_[i],2.0),0.5);
    if(i < bpm_.GetNumberOfSteps() / 2) { // first half
      if(horizontal_first){
        displacement_.push_back(std::make_pair(displacement,1));
      } else {
        displacement_.push_back(std::make_pair(displacement,0));
      }
    } else { // second half
      if(horizontal_first){
        displacement_.push_back(std::make_pair(displacement,0));
      } else {
        displacement_.push_back(std::make_pair(displacement,1));
      }
    }
  }
  int step_i = 1;
  for( auto i = displacement_.begin(); i != displacement_.end(); ++i) {
    auto displacement = (*i).first;
    auto scan = (*i).second;
    std::string scan_name = "";
    if( scan ) { 
      scan_name = "Horizontal";
    } else {
      scan_name = "Vertical";
    }
    std::cout << "Step: " << step_i << " , Displacement: " << displacement << ", Scan: " << scan_name << std::endl;
    step_i++;
  }
  LoadSteps();
  InitHistograms();
  return 0;
}

int HourglassData::LoadSteps() {
  std::ifstream in_file(epoch_steps_file_name_.c_str());
  std::string line;
  if(in_file) {
    while(getline(in_file,line)) {
      if(line[0] == '#') continue;
      std::stringstream ss;
      ss.str(line);
      long int start_time;
      long int end_time;
      ss >> start_time >> end_time;
      step_boundaries_.push_back(std::make_pair(start_time,end_time));
    }
  } else {
    std::cerr << "could not open " << epoch_steps_file_name_ << std::endl;
    return 1;
  } 
  return 0;
}

int HourglassData::InitHistograms() {
  int nbins = 100;
  float min_zvtx = -300;
  float max_zvtx = 300;
  for(unsigned int i = 0; i < step_boundaries_.size(); i++){
    auto displacement = displacement_[i].first;
    auto scan         = displacement_[i].second;
    std::string scan_name = "";
    if( scan ) { 
      scan_name = "Horizontal";
    } else {
      scan_name = "Vertical";
    }
    std::stringstream name;
    std::stringstream title;
    name << "bbc_zvtx_step_" << i;
    title << scan_name << " scan, " << displacement << " microns step, Run : " << run_number_ << ", BBC Z-vertex distribution;z-vertex;counts";
    TH1F* bbc = new TH1F(name.str().c_str(),title.str().c_str(),nbins,min_zvtx,max_zvtx);
    bbc->Sumw2();

    name.str("");
    title.str("");
    name << "zdc_zvtx_step_" << i;
    title << scan_name << " scan, " << displacement << " microns step, Run : " << run_number_ << ", ZDC Z-vertex distribution;z-vertex;counts";
    TH1F* zdc = new TH1F(name.str().c_str(),title.str().c_str(),nbins,min_zvtx,max_zvtx);
    zdc->Sumw2();

    bbc_z_vtx_.insert( std::pair<int,TH1F*>(i,bbc) );
    zdc_z_vtx_.insert( std::pair<int,TH1F*>(i,zdc) );
    plot_registry_.push_back(bbc);
    plot_registry_.push_back(zdc);
  }
  return 0;
}

int HourglassData::Run() {
  // Run 12 Defaults
  int bbc_wide_trig_bit = 0x00000002;
  int zdc_wide_trig_bit = 0x00000004;

  std::ifstream in_file(scalers_file_name_.c_str());
  std::string line;
  if(in_file) {
    while(getline(in_file,line)) {
      if(line[0] == '#') continue;
      std::stringstream ss;
      ss.str(line);
      long int event_number;
      MergedPRDFDST d;
      d.Reset();

      ss >> event_number 
        >> d.atp_number 
        >> d.timestamp 
        >> d.bunch 
        >> d.gl1p_bbc 
        >> d.gl1p_clock 
        >> d.gl1p_zdc_wide 
        >> d.gl1p_zdc_narrow 
        >> d.trigraw    
        >> d.triglive   
        >> d.trigscaled 
        >> d.bbc_z      
        >> d.zdc_z;

      int step = bpm_.LookupStep((double)d.timestamp);
        auto bbc = bbc_z_vtx_.find(step);
        auto zdc = zdc_z_vtx_.find(step);
        if( (bbc != bbc_z_vtx_.end() ) && (zdc != zdc_z_vtx_.end()) ) { 
          if( ( d.trigscaled & bbc_wide_trig_bit ) > 0) { (bbc->second)->Fill(d.bbc_z); }
          if( ( d.trigscaled & zdc_wide_trig_bit ) > 0) { (zdc->second)->Fill(d.zdc_z); }
        }
    }
  } else {
    std::cerr << "could not open " << scalers_file_name_ << std::endl;
    return 1;
  } 
  return 0;
}

int HourglassData::CreatePlots() {

 return 0;
}

int HourglassData::SaveFigures( const std::string& figure_output_dir = "./") {
  std::string tfile_name = figure_output_dir + "/" + run_number_ + "_HourglassData.root"; 
  std::string name       = figure_output_dir + "/" + run_number_ + "_HourglassData.pdf";
  std::string out_file_name = name + "[";
  TCanvas* booklet = new TCanvas("booklet","BBC Efficiency Plots");
  TFile* root_out = new TFile(tfile_name.c_str(), "RECREATE");
  booklet->Print(out_file_name.c_str());
  for(auto plot_i = plot_registry_.begin(); plot_i != plot_registry_.end(); ++plot_i) {
    auto draw_obj = *plot_i;
    if(!draw_obj) continue;
    std::cout << "Saving/Drawing: " << draw_obj->GetName() << std::endl;
    root_out->cd();
    draw_obj->Write();
    if(!draw_obj) continue;
    TCanvas* c = new TCanvas();
    c->cd();
    draw_obj->Draw();
    c->Print(name.c_str());
    delete c;
  }
  out_file_name = name + "]";
  booklet->Print(out_file_name.c_str());
  root_out->Write();
  root_out->Close();
  if(root_out) delete root_out;
  delete booklet;
  return 0;
}
