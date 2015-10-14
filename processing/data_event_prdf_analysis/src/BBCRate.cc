#include "BBCRate.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <string>
#include <iomanip>
#include <utility>
#include <cmath>

#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"

const double BBCRate::CLOCK_RATE = 9.5e6;

BBCRate::BBCRate() {
  std::cout << "Creating BBCRate instance at: " << this << std::endl;
}

BBCRate::~BBCRate() {
  std::cout << "Destroying instance of BBCRate at: " << this << std::endl;
  for(auto i = plot_registry_.begin(); i != plot_registry_.end(); ++i) {
    auto deleteme = *i;
    std::string name = deleteme->GetName();
    if(deleteme) {
      delete deleteme;
    //  std::cout << "Cleaning up: " << name << std::endl;
    }
  }
}


int BBCRate::Init(const std::string& gl1p_scalers_data_file, const std::string& run_number_ ) {  
  std::ifstream in_file(gl1p_scalers_data_file.c_str());            
  this->run_number_ = run_number_;
  // get min/max just in case we want to draw a frame later on
  bool first_line = true;
  if( in_file ) {
    std::string line = "";
    while(getline(in_file,line)) {
      if(line[0] == '#') continue;
      //std::cout << line << std::endl;
      long int time = -999;
      long int gl1p_scaler_bbc = -999;
      long int gl1p_scaler_clock = -999;
      long int gl1p_scaler_zdc_wide = -999;
      long int gl1p_scaler_zdc_narrow = -999;
      // get min and max rates for frame-setting or other drawing needs
      if(first_line){ 
        bbc_max_rate_ = GetGL1PRate(gl1p_scaler_bbc,gl1p_scaler_clock);
        bbc_min_rate_ = GetGL1PRate(gl1p_scaler_bbc,gl1p_scaler_clock);
        first_line = false;
      } else {
        double test_rate = GetGL1PRate(gl1p_scaler_bbc,gl1p_scaler_clock);
        if (test_rate > bbc_max_rate_ ) {bbc_max_rate_ = test_rate;}
        if (test_rate < bbc_min_rate_ ) {bbc_min_rate_ = test_rate;}
      }
      std::stringstream ss;
      ss << line;
      ss >> time 
       >> gl1p_scaler_bbc 
       >> gl1p_scaler_clock 
       >> gl1p_scaler_zdc_wide 
       >> gl1p_scaler_zdc_narrow ;
      GL1P gl1p;
      gl1p.Init();
      gl1p.bbc_ = gl1p_scaler_bbc;
      gl1p.clock_ = gl1p_scaler_clock;
      gl1p.zdc_w_ = gl1p_scaler_zdc_wide;
      gl1p.zdc_n_ = gl1p_scaler_zdc_narrow;
      // std::cout << time << " " << gl1p.GetScalerSerializationString() << std::endl;
      gl1p_scalers_[time] = gl1p;
    }
  } else {
    std::cerr << gl1p_scalers_data_file << " did not open correctly." << std::endl;
    return 1;
  }
  std::cout << "Loaded: " << gl1p_scalers_.size() 
      << " data points from: " << gl1p_scalers_data_file << std::endl;
  time_offset_ = gl1p_scalers_.begin()->first;
  return 0;
}

int BBCRate::MakeHistograms() {
  long int low_range = gl1p_scalers_.begin()->first - time_offset_;
  long int high_range = (--gl1p_scalers_.end())->first - time_offset_;
  long int bins = high_range - low_range;
  double low = (double)low_range;
  double high = (double)high_range;
  low -=0.5;
  high-=0.5;
  std::stringstream title;

  // Remember, for scaler data, live time is an issue. Ratios of clock to relevant scalar
  // will produce the appropriate rate, if you correct for clock tick rate.
  title << "BBC Rate, Run: " << run_number_ << "; time(s); BBCLL1(>0 tubes) GL1P Scaler";
  bbc_rate_ = new TH1F("bbc_rate",title.str().c_str(),bins,low,high); 
  plot_registry_.push_back(bbc_rate_);
  title.str("");
  title << "ZCD Narrow Rate, Run: " << run_number_ << "; time(s); ZDCLL1Narrow GL1P Scaler";
  zdc_n_rate_ = new TH1F("zdc_n_rate",title.str().c_str(),bins, low, high);
  plot_registry_.push_back(zdc_n_rate_);
  title.str("");
  title << "ZDC Wide Rate, Run: " << run_number_ << "; time(s); ZDCLL1Wide GL1P Scaler";
  zdc_w_rate_ = new TH1F("zdc_w_rate",title.str().c_str(),bins, low, high);
  plot_registry_.push_back(zdc_w_rate_);
  title.str("");
  title << "Clock Rate, Run: " << run_number_ << "; time(s); Clock GL1P Scaler";
  clock_rate_ = new TH1F("clock_rate", title.str().c_str(), bins, low, high);
  plot_registry_.push_back(clock_rate_);
  title.str("");

  for(auto i = gl1p_scalers_.begin(); i != gl1p_scalers_.end(); ++i) {
    double time = (double) i->first;
    time = time - time_offset_;
    GL1P gl1p = i->second;
    if(gl1p.clock_ == 0) {
      std::cout << "Zero clock counts, can't calculate rate at time: " << time << std::endl;
      continue;
    }
    
    int lookup_bin = bbc_rate_->FindBin(time);
    bbc_rate_data_[i->first];
    bbc_rate_data_[i->first].Reset();
    bbc_rate_data_[i->first].rate  = GetGL1PRate(gl1p.bbc_,gl1p.clock_);
    bbc_rate_data_[i->first].rate_err = GetGL1PRateError(gl1p.bbc_,gl1p.clock_);
    bbc_rate_data_[i->first].bbc_gl1p = gl1p.bbc_;
    bbc_rate_data_[i->first].clock_gl1p = gl1p.clock_;

    bbc_rate_->SetBinContent(lookup_bin, GetGL1PRate(gl1p.bbc_,gl1p.clock_));
    bbc_rate_->SetBinError  (lookup_bin, GetGL1PRateError(gl1p.bbc_,gl1p.clock_));
    zdc_w_rate_->SetBinContent(lookup_bin, GetGL1PRate(gl1p.zdc_w_,gl1p.clock_));
    zdc_w_rate_->SetBinError  (lookup_bin, GetGL1PRateError(gl1p.zdc_w_,gl1p.clock_));
    zdc_n_rate_->SetBinContent(lookup_bin, GetGL1PRate(gl1p.zdc_n_,gl1p.clock_));
    zdc_n_rate_->SetBinError  (lookup_bin, GetGL1PRateError(gl1p.zdc_n_,gl1p.clock_));
    clock_rate_->SetBinContent(lookup_bin, (double)gl1p.clock_);
    clock_rate_->SetBinError(lookup_bin, pow((double)gl1p.clock_,0.5));
  }
  return 0;
}

double BBCRate::GetGL1PRate(long int gl1p_scaler, long int gl1p_clock){
  double rate = 0.;
  double n1 = (double)gl1p_scaler;
  double n2 = (double)gl1p_clock;
  if( gl1p_clock == 0 ){
    std::cout << "Rates are meaningless when divided by zero... clock scaler is zero..."
        << std::endl
        << " rate is set to zero, but you shouldn't be seeing this error with proper use..."
        << std::endl;
    return 0.;
  } else if (gl1p_scaler == 0 ) {
    return 0.0;
  } else {
    rate = n1/n2*CLOCK_RATE; 
  }
  return rate;
}

double BBCRate::GetGL1PRateError(long int gl1p_scaler, long int gl1p_clock){
  double error = 0.;
  double rate = GetGL1PRate(gl1p_scaler,gl1p_clock);
  double n1 = (double)gl1p_scaler;
  double n2 = (double)gl1p_clock;
  double n1_err = pow((double)gl1p_scaler,0.5);
  double n2_err = pow((double)gl1p_clock ,0.5);
  if( gl1p_scaler == 0 ) {
    return 0.;
  } else if (gl1p_clock == 0 ) {
    std::cout << "Rates are meaningless when divided by zero... clock scaler is zero..."
        << std::endl
        << " rate is set to zero, but you shouldn't be seeing this error with proper use..."
        << std::endl;
    return 0.;
  } else {
    error = rate*pow( pow(n1_err/n1,2.) + pow(n2_err/n2,2.) , 0.5 );
  }
  return error;
}

int BBCRate::Run() {
  MakeHistograms();
  return 0;
}

int BBCRate::MakeFigures(const std::string& figure_output_dir) {
  std::string tfile_name = figure_output_dir + "/" + run_number_ + "_BBCRatePlots.root"; 
  std::string pdf_name   = figure_output_dir + "/" + run_number_ + "_BBCRatePlots.pdf";
  std::string out_file_name = pdf_name + "[";
  TCanvas* booklet = new TCanvas("booklet","Scaler Rate Plots");
  TFile* root_out = new TFile(tfile_name.c_str(), "RECREATE");
  booklet->Print(out_file_name.c_str());
  for(auto plot_i = plot_registry_.begin(); plot_i != plot_registry_.end(); ++plot_i) {
    auto draw_obj = *plot_i;
    std::string name = draw_obj->GetName();
    if(!draw_obj) continue;
    //std::cout << "Saving/Drawing: " << name << std::endl;
    root_out->cd();
    draw_obj->Write();
    if(!draw_obj) continue;
    TCanvas* c = new TCanvas();
    c->cd();
    draw_obj->Draw();
    c->Print(pdf_name.c_str()); // ->GetName() returns const char*
    delete c;
  }
  out_file_name = pdf_name + "]";
  booklet->Print(out_file_name.c_str());
  root_out->Write();
  root_out->Close();
  delete booklet;
  return 0;
}
