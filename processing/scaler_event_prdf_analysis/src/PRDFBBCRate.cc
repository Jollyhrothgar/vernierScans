#include "PRDFBBCRate.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <string>
#include <iomanip>
#include <utility>

#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"


PRDFBBCRate::PRDFBBCRate() {
  std::cout << "Creating PRDFBBCRate instance at: " << this << std::endl;
}

PRDFBBCRate::~PRDFBBCRate() {
  std::cout << "Destroying instance of PRDFBBCRate at: " << this << std::endl;
}


int PRDFBBCRate::Init(
  const std::string& run_number, 
  const std::string& prdf_data_file
) {  
  run_number_ = run_number;
  bbc_raw_rate_vs_time_ = new TGraph();
  plot_registry_.push_back(bbc_raw_rate_vs_time_);
  bbc_raw_rate_vs_time_ -> SetName("bbc_raw_rate_vs_time_");
  bbc_raw_rate_vs_time_ -> SetTitle("Raw BBC Rate vs Time From Scaler Events; Time (s); BBC Rate (Hz)");
  bbc_raw_rate_vs_time_ -> SetMarkerStyle(kFullCircle);
  bbc_raw_rate_vs_time_ -> SetMarkerSize(0.75);
  bbc_raw_rate_vs_time_ -> SetMarkerColor(kRed);
  bbc_raw_rate_vs_time_ -> SetLineColor(kBlack);

  scaler_vs_time_["bbc_w_raw"] = new TGraph();
  scaler_vs_time_["bbc_n_raw"] = new TGraph();
  scaler_vs_time_["bbc_w_live"] = new TGraph();
  scaler_vs_time_["bbc_n_live"] = new TGraph();
  scaler_vs_time_["clock_raw"] = new TGraph();
  scaler_vs_time_["clock_live"] = new TGraph();
  scaler_vs_time_["zdc_w_raw"] = new TGraph();
  scaler_vs_time_["zdc_w_live"] = new TGraph();
  scaler_vs_time_["live_time_bbc_w"] = new TGraph();
  scaler_vs_time_["live_time_bbc_n"] = new TGraph();
  scaler_vs_time_["live_time_clock"] = new TGraph();
  scaler_vs_time_["live_time_zdc_w"] = new TGraph();
  scaler_vs_time_["i_live_time_bbc_w"] = new TGraph();
  scaler_vs_time_["i_live_time_bbc_n"] = new TGraph();
  scaler_vs_time_["i_live_time_clock"] = new TGraph();
  scaler_vs_time_["i_live_time_zdc_w"] = new TGraph();
  
  plot_registry_.push_back(scaler_vs_time_["bbc_w_raw"]);
  plot_registry_.push_back(scaler_vs_time_["bbc_n_raw"]);
  plot_registry_.push_back(scaler_vs_time_["bbc_w_live"]);
  plot_registry_.push_back(scaler_vs_time_["bbc_n_live"]);
  plot_registry_.push_back(scaler_vs_time_["clock_raw"]);
  plot_registry_.push_back(scaler_vs_time_["clock_live"]);
  plot_registry_.push_back(scaler_vs_time_["zdc_w_raw"]);
  plot_registry_.push_back(scaler_vs_time_["zdc_w_live"]);
  plot_registry_.push_back(scaler_vs_time_["live_time_bbc_w"]);
  plot_registry_.push_back(scaler_vs_time_["live_time_bbc_n"]);
  plot_registry_.push_back(scaler_vs_time_["live_time_clock"]);
  plot_registry_.push_back(scaler_vs_time_["live_time_zdc_w"]);
  plot_registry_.push_back(scaler_vs_time_["i_live_time_bbc_w"]);
  plot_registry_.push_back(scaler_vs_time_["i_live_time_bbc_n"]);
  plot_registry_.push_back(scaler_vs_time_["i_live_time_clock"]);
  plot_registry_.push_back(scaler_vs_time_["i_live_time_zdc_w"]);

  for(auto gr_itr = scaler_vs_time_.begin(); gr_itr != scaler_vs_time_.end(); ++gr_itr) {
    auto graph = gr_itr->second;
    std::string name = gr_itr->first;
    std::string title = run_number+ ": " + name + " vs time; seconds;" + name;

    graph->SetName(name.c_str());
    graph->SetTitle(title.c_str());
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);
  }
  std::ifstream in_file(prdf_data_file.c_str());            
  if( in_file ) {
    std::string line;
    while(getline(in_file,line)) {
      if(line[0] == '#') continue;
      //std::cout << line << std::endl;
      std::stringstream ss;
      ss.str(line.c_str());
      unsigned long long epoch;
      Scaler s;
      s.Init();
      ss >> epoch 
          >> s.clock_raw_  >> s.bbc_w_raw_  >> s.bbc_n_raw_  >> s.zdc_w_raw_
          >> s.clock_live_ >> s.bbc_w_live_ >> s.bbc_n_live_ >> s.zdc_w_live_;
      //std::cout << epoch << " " << s.GetScalerSerializationString() << std::endl;
      prdf_data_[epoch] = s;
    }
  } else {
    std::cerr << prdf_data_file << " did not open correctly." << std::endl;
    return 1;
  }
  std::cout << "Loaded: " << prdf_data_.size() 
      << " data points from: " << prdf_data_file << std::endl;

  /* Track clock ticks per second..I think
  long double num = second->second.clock_raw_;
  long double denom = (long double)(second->first - first->first);
  */
  auto final = prdf_data_.end();
  final--;
  std::cout << "Raw Trigger   \"BBCLL1(>0 tubes) novertex\" : " << final->second.bbc_w_raw_  << std::endl;
  std::cout << "Raw Trigger   \"BBCLL1(>0 tubes)\"          : " << final->second.bbc_n_raw_  << std::endl;
  std::cout << "Raw Trigger   \"ZDCLL1wide\"                : " << final->second.zdc_w_raw_  << std::endl;
  std::cout << "Raw Trigger   \"Clock\"                    : " << final->second.clock_raw_  << std::endl;
  std::cout << "Live Trigger: \"BBCLL1(>0 tubes) novertex\" : " << final->second.bbc_w_live_ << std::endl;
  std::cout << "Live Trigger: \"BBCLL1(>0 tubes)\"          : " << final->second.bbc_n_live_ << std::endl;
  std::cout << "Live Trigger: \"ZDCLL1wide\"                : " << final->second.zdc_w_live_ << std::endl;
  std::cout << "Live Trigger: \"Clock\"                     : " << final->second.clock_live_ << std::endl;
  return 0;
}

int PRDFBBCRate::Run() {
  auto start_itr = prdf_data_.begin();
  start_itr++;
  for(auto evt_i = start_itr; evt_i != prdf_data_.end(); ++evt_i) {
    //std::cout << evt_i->first - (prdf_data_.begin()->first) << std::endl;
    auto next = evt_i;
    ++next;
    if( next != prdf_data_.end() ) {
      long double elapsed_time      = (double)(next->first - evt_i->first);
      long double bbc_w_raw_counts  = (double)next->second.bbc_w_raw_   - (double)evt_i->second.bbc_w_raw_;
      long double bbc_w_live_counts = (double)next->second.bbc_w_live_  - (double)evt_i->second.bbc_w_live_;
      long double bbc_n_raw_counts  = (double)next->second.bbc_n_raw_   - (double)evt_i->second.bbc_n_raw_;
      long double bbc_n_live_counts = (double)next->second.bbc_n_live_  - (double)evt_i->second.bbc_n_live_;
      long double clock_raw_counts  = (double)next->second.clock_raw_   - (double)evt_i->second.clock_raw_;
      long double clock_live_counts = (double)next->second.clock_live_  - (double)evt_i->second.clock_live_;
      long double zdc_w_raw_counts  = (double)next->second.zdc_w_raw_   - (double)evt_i->second.zdc_w_raw_;
      long double zdc_w_live_counts = (double)next->second.zdc_w_live_  - (double)evt_i->second.zdc_w_live_;

      //std::cout << bbc_w_raw_counts << ", " << elapsed_time << std::endl;
      double approx_time = (double)(evt_i->first - prdf_data_.begin()->first) + elapsed_time/2.0;
      double real_time = (double)(evt_i->first) + elapsed_time/2.0;
      // double bbc_rate = bbc_w_raw_counts / elapsed_time;
      double test_rate = bbc_n_live_counts/clock_live_counts*9.5e6;

      // std::cout << "Live rate, BBCLL1(>0 tubes): " << std::setw(12) << bbc_n_live_counts / elapsed_time  << " Hz " << std::endl;
      // std::cout << "Test rate, BBCLL1(>0 tubes): " << std::setw(12) << test_rate << " Hz " << std::endl;
      float bbc_w_livetime = bbc_w_live_counts/bbc_w_raw_counts;

      // Already have a graph called g_bbc_w_live
      g_bbc_w_live->SetPoint(g_bbc_w_live->GetN(), time, bbc_w_livetime);

      float bbc_n_livetime = bbc_n_live_counts/bbc_n_raw_counts;
      float clock_livetime = clock_live_counts/clock_raw_counts;
      float zdc_w_livetime = zdc_w_live_counts/zdc_w_raw_counts;

      livetime_[real_time].bbc_w = bbc_w_livetime;
      livetime_[real_time].bbc_n = bbc_n_livetime;
      livetime_[real_time].clock = clock_livetime;
      livetime_[real_time].zdc_w = zdc_w_livetime;

      std::cout << std::setprecision(12) << real_time << " " << livetime_[real_time].GetScalerSerializationString() << std::endl;
      
      bbc_raw_rate_vs_time_->SetPoint(bbc_raw_rate_vs_time_->GetN(),approx_time, test_rate);

      bbc_w_raw_rate_.push_back(std::make_pair(real_time,bbc_w_raw_counts/elapsed_time)); 
      bbc_w_live_rate_.push_back(std::make_pair(real_time,bbc_w_live_counts/elapsed_time));
      bbc_n_raw_rate_.push_back(std::make_pair(real_time,bbc_n_raw_counts/elapsed_time));
      bbc_n_live_rate_.push_back(std::make_pair(real_time,bbc_n_live_counts/elapsed_time));
      scaler_vs_time_["bbc_w_raw" ] ->SetPoint( scaler_vs_time_["bbc_w_raw" ] ->GetN(), approx_time, bbc_w_raw_counts  /elapsed_time);
      scaler_vs_time_["bbc_w_live"] ->SetPoint( scaler_vs_time_["bbc_w_live"] ->GetN(), approx_time, bbc_w_live_counts /elapsed_time);
      scaler_vs_time_["bbc_n_raw" ] ->SetPoint( scaler_vs_time_["bbc_n_raw" ] ->GetN(), approx_time, bbc_n_raw_counts  /elapsed_time);
      scaler_vs_time_["bbc_n_live"] ->SetPoint( scaler_vs_time_["bbc_n_live"] ->GetN(), approx_time, bbc_n_live_counts /elapsed_time);
      scaler_vs_time_["clock_raw" ] ->SetPoint( scaler_vs_time_["clock_raw" ] ->GetN(), approx_time, clock_raw_counts  /elapsed_time);
      scaler_vs_time_["clock_live"] ->SetPoint( scaler_vs_time_["clock_live"] ->GetN(), approx_time, clock_live_counts /elapsed_time);
      scaler_vs_time_["zdc_w_raw" ] ->SetPoint( scaler_vs_time_["zdc_w_raw" ] ->GetN(), approx_time, zdc_w_raw_counts  /elapsed_time);
      scaler_vs_time_["zdc_w_live"] ->SetPoint( scaler_vs_time_["zdc_w_live"] ->GetN(), approx_time, zdc_w_live_counts /elapsed_time);
  
      scaler_vs_time_["live_time_bbc_w"]->SetPoint( scaler_vs_time_["live_time_bbc_w"]->GetN(), approx_time, bbc_w_livetime );
      scaler_vs_time_["live_time_bbc_n"]->SetPoint( scaler_vs_time_["live_time_bbc_n"]->GetN(), approx_time, bbc_n_livetime );
      scaler_vs_time_["live_time_clock"]->SetPoint( scaler_vs_time_["live_time_clock"]->GetN(), approx_time, clock_livetime );
      scaler_vs_time_["live_time_zdc_w"]->SetPoint( scaler_vs_time_["live_time_zdc_w"]->GetN(), approx_time, zdc_w_livetime );
    }
  }
  double first_time = livetime_.begin()->first;
  auto end = prdf_data_.end();
  --end;
  double last_time = end->first;
  for(double time = first_time; time < last_time; time+=0.5) {
    double left_time;
    double right_time;
    // float bbc_w_livetime;
    // float bbc_n_livetime;
    // float clock_livetime;
    // float zdc_w_livetime;
    if(livetime_.find(time) != livetime_.end()) { // A real live-time value exists
      std::cout << time << std::endl;
      scaler_vs_time_["i_live_time_bbc_w"]->SetPoint(scaler_vs_time_["i_live_time_bbc_w"]->GetN(), time - first_time, livetime_[time].bbc_w );
      scaler_vs_time_["i_live_time_bbc_n"]->SetPoint(scaler_vs_time_["i_live_time_bbc_n"]->GetN(), time - first_time, livetime_[time].bbc_n );
      scaler_vs_time_["i_live_time_clock"]->SetPoint(scaler_vs_time_["i_live_time_clock"]->GetN(), time - first_time, livetime_[time].clock );
      scaler_vs_time_["i_live_time_zdc_w"]->SetPoint(scaler_vs_time_["i_live_time_zdc_w"]->GetN(), time - first_time, livetime_[time].zdc_w );
    } else { // we must intripolate the live-time 
     FindLivetimeBounds(time,left_time,right_time);
     float bbc_w_livetime_left  = livetime_[left_time].bbc_w;
     float bbc_n_livetime_left  = livetime_[left_time].bbc_n;
     float clock_livetime_left  = livetime_[left_time].clock;
     float zdc_w_livetime_left  = livetime_[left_time].zdc_w;
     float bbc_w_livetime_right = livetime_[right_time].bbc_w;
     float bbc_n_livetime_right = livetime_[right_time].bbc_n;
     float clock_livetime_right = livetime_[right_time].clock;
     float zdc_w_livetime_right = livetime_[right_time].zdc_w;

     std::cout << "bbc_w_livetime_left   " << bbc_w_livetime_left  << std::endl;
     std::cout << "bbc_n_livetime_left   " << bbc_n_livetime_left  << std::endl;
     std::cout << "clock_livetime_left   " << clock_livetime_left  << std::endl;
     std::cout << "zdc_w_livetime_left   " << zdc_w_livetime_left  << std::endl;
     std::cout << "bbc_w_livetime_right  " << bbc_w_livetime_right << std::endl;
     std::cout << "bbc_n_livetime_right  " << bbc_n_livetime_right << std::endl;
     std::cout << "clock_livetime_right  " << clock_livetime_right << std::endl;
     std::cout << "zdc_w_livetime_right  " << zdc_w_livetime_right << std::endl;
    // scaler_vs_time_["i_live_time_bbc_w"]->SetPoint(scaler_vs_time_["i_live_time_bbc_w"]->GetN(), time - first_time, IntripolateLivetime(time,left_time,bbc_w_livetime_left,right_time,bbc_w_livetime_right));
    // scaler_vs_time_["i_live_time_bbc_n"]->SetPoint(scaler_vs_time_["i_live_time_bbc_n"]->GetN(), time - first_time, IntripolateLivetime(time,left_time,bbc_n_livetime_left,right_time,bbc_n_livetime_right));
    // scaler_vs_time_["i_live_time_clock"]->SetPoint(scaler_vs_time_["i_live_time_clock"]->GetN(), time - first_time, IntripolateLivetime(time,left_time,clock_livetime_left,right_time,clock_livetime_right));
    // scaler_vs_time_["i_live_time_zdc_w"]->SetPoint(scaler_vs_time_["i_live_time_zdc_w"]->GetN(), time - first_time, IntripolateLivetime(time,left_time,zdc_w_livetime_left,right_time,zdc_w_livetime_right));
    }
  }
  return 0;
}
bool PRDFBBCRate::FindLivetimeBounds(
    double target_time,
    double& left_time,
    double& right_time 
    ){
  auto low = prdf_data_.lower_bound(target_time);
  if (low == prdf_data_.end()) {
    // debug wp, this may need work, because it will lead to using non-existent keys.
    // std::cout << "Lower bound lookup not found, skip this point." << std::endl;
    return false;
  } else if (low == prdf_data_.begin()) {
    std::cout << "Found lowest point" << std::endl;
    return true;
  } else {
    auto prev = low;
    --prev;
    if ((target_time - prev->first) < (low->first - target_time)) {
      // debug
      // std::cout << "lower bound = " << prev->first << std::endl;
      // std::cout << prev->first - first_time_;
      // prev++;
      // std::cout << ", " << target_time - first_time_ << ", " << prev->first - first_time_ << std::endl;
      left_time = prev->first;
      prev++;
      right_time = prev->first;
      return true;
    } else {
      // debug
      // std::cout << "lower bound = " << low->first << std::endl;
      // std::cout << prev->first - first_time_;
      // prev++;
      // std::cout << ", " << target_time - first_time_ << ", " << prev->first - first_time_ << std::endl;
      left_time = prev->first;
      prev++;
      right_time = prev->first;
      return true;
    }
  }
  return true;
}

float PRDFBBCRate::IntripolateLivetime(
  double target_time,
  double left_time,
  double left_livetime,
  double right_time,
  double right_livetime
){
  double lt = left_time;
  double rt = right_time;
  double ll = left_livetime;
  double rl = right_livetime;
  double tt = target_time;

  long double slope = (rl - ll)/(rt - lt);

  long double diff = tt - lt;
  return diff*slope + rt;
}

int PRDFBBCRate::SaveFigures(const std::string& figure_output_dir) {
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
