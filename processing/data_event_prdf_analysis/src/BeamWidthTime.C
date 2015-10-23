#include "BeamWidthTime.h"

#include "TGraph.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TError.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

const double BeamWidthTime::CLOCK_RATE = 9.5e6;

BeamWidthTime::BeamWidthTime(){
  std::stringstream ss;
  ss << "BeamWidthTime" << "_0x" << std::hex << this << std::endl;
  this_name_ = ss.str();
  std::cout << "Creating instance: " << this_name_ << std::endl;
  global_output_.str("");
}

BeamWidthTime::~BeamWidthTime() {
  std::cout << "Destroying instance: " << this_name_ << std::endl;
  std::cout << global_output_.str() << std::endl;
  for(auto i = plot_registry_.begin(); i != plot_registry_.end(); ++i) {
    auto deleteme = *i;
    std::string name = deleteme->GetName();
    if(deleteme) {
      delete deleteme;
      // std::cout << "Cleaning up: " << name << std::endl;
    }
  }
}

int BeamWidthTime::Init(
    const std::string& run_number ,
    const std::string& bpm_data   ,
    const std::string& bpm_steps  ,
    const std::string& epoch_steps,
    const std::string& scalers_file,
    const std::string& planned_beam_steps_file
    ){
  bbc.Init(
      scalers_file,
      run_number);
  bbc.Run();
  bpm.Init(
      run_number,
      bpm_data,
      bpm_steps
    );
  bpm.Run();

  run_number_ = run_number;
  bbc_data_ = bbc.GetBBCRateData();
  bpm_data_ = bpm.GetBeamSeparationData();

  step_time_ = new TH1F("step_time_","Average Time Per Step",20,0,120); 
  plot_registry_.push_back(step_time_);

  // load boundaries
  std::ifstream in_file(epoch_steps.c_str());
  if(in_file) {
    std::string line;
    while(getline(in_file,line)) {
      std::stringstream ss;
      long int start;
      long int end;
      ss << line;
      if( line[0] == '#' ) continue;
      ss >> start >> end;
      step_time_->Fill((double)end - (double)start);
      steps_epoch_boundaries_.push_back(std::make_pair(start,end));
    }
  } else {
    std::cout << "BeamWidthTime did could not load epoch step boundaries. Segfault comming soon! " << std::endl;
    return 1;
  }
  global_output_ << "There were " << steps_epoch_boundaries_.size() << " scanned steps." << std::endl;
  global_output_ << "  Time per step: " << step_time_->GetMean() << " seconds." << std::endl;

  /** Generate Scan Boundaries */
  long int first_scan_start = steps_epoch_boundaries_[0].first;
  long int first_scan_end   = steps_epoch_boundaries_[(steps_epoch_boundaries_.size())/2-1].second;
  long int second_scan_start = steps_epoch_boundaries_[(steps_epoch_boundaries_.size())/2].first;
  long int second_scan_end   = steps_epoch_boundaries_[steps_epoch_boundaries_.size() - 1].second;
  scan_1_ = std::make_pair( first_scan_start  , first_scan_end  );
  scan_2_ = std::make_pair( second_scan_start , second_scan_end );

  /** Initialize the rest of the plots */ 
  std::stringstream name;
  std::stringstream title;
  rate_vs_sep_s0_h = new TGraphErrors();
  rate_vs_sep_s0_h->SetName("scan0_sep0");
  rate_vs_sep_s0_v = new TGraphErrors();
  rate_vs_sep_s0_v->SetName("scan0_sep1");
  rate_vs_sep_s1_h = new TGraphErrors();
  rate_vs_sep_s1_h->SetName("scan1_sep0");
  rate_vs_sep_s1_v = new TGraphErrors();
  rate_vs_sep_s1_v->SetName("scan1_sep1");
  rate_vs_sep_s0_h->SetTitle("Separation vs BBC Rate, Scan 0, Horizontal;Separation (microns); BBC Rate (Hz)");
  rate_vs_sep_s0_v->SetTitle("Separation vs BBC Rate, Scan 0, Vertical  ;Separation (microns); BBC Rate (Hz)");
  rate_vs_sep_s1_h->SetTitle("Separation vs BBC Rate, Scan 1, Horizontal;Separation (microns); BBC Rate (Hz)");
  rate_vs_sep_s1_v->SetTitle("Separation vs BBC Rate, Scan 1, Vertical  ;Separation (microns); BBC Rate (Hz)");
  rate_vs_sep_s0_h->SetMarkerStyle(7);
  rate_vs_sep_s0_v->SetMarkerStyle(7);
  rate_vs_sep_s1_h->SetMarkerStyle(7);
  rate_vs_sep_s1_v->SetMarkerStyle(7);
  plot_registry_.push_back(rate_vs_sep_s0_h);
  plot_registry_.push_back(rate_vs_sep_s0_v);
  plot_registry_.push_back(rate_vs_sep_s1_h);
  plot_registry_.push_back(rate_vs_sep_s1_v);

  /** Load Planned Beam Steps File */
  std::ifstream in_planned_steps(planned_beam_steps_file.c_str());
  if(in_planned_steps) {
    std::string line;
    while(getline(in_planned_steps,line)) {
      if(line[0] == '#') continue; 
      std::stringstream ss;
      float step;
      ss << line;
      ss >> step;
      planned_beam_separation_.push_back(step*1000.0); // file contains steps in millimeters, convert to micrometers
    }
  } else {
    std::cout << "Coule not properly load the planned beam steps from file: " << planned_beam_steps_file << std::endl;
    return 1;
  }
  std::cout << "Got the " << planned_beam_separation_.size() << " planned steps from " << planned_beam_steps_file << std::endl;
  return 0;
}

int BeamWidthTime::CorrelateTime() {
  BBCRateData bbc;
  bbc.Reset();
  BeamSeparationData bpm;
  bpm.Reset();
  double sum_h1 = 0.;
  double sum_h2 = 0.;
  for ( auto i = bbc_data_.begin(); i != bbc_data_.end(); ++i) { // Lookup loop
    bbc.Reset();
    bpm.Reset();
    BBCRateData bbc = i->second;
    long int time = i->first;
    if(!IsStep(time)) continue;
    auto bpm_lookup = bpm_data_.find(time);
    if( bpm_lookup != bpm_data_.end() ) { // found time index!
      BeamSeparationData bpm = bpm_lookup->second;
      if( time >= scan_1_.first && time <= scan_1_.second ) {
        // first scan
        TGraphErrors* h  = rate_vs_sep_s0_h;
        TGraphErrors* v  = rate_vs_sep_s0_v;
        h->SetPoint(h->GetN(), bpm.x, bbc.rate);
        h->SetPointError(h->GetN()-1, bpm.x_err, bbc.rate_err);
        v->SetPoint(v->GetN(), bpm.y, bbc.rate);
        v->SetPointError(v->GetN()-1, bpm.y_err, bbc.rate_err);
        sum_h1 += fabs(bpm.x);
      } else { 
        // second scan
        TGraphErrors* h  = rate_vs_sep_s1_h;
        TGraphErrors* v  = rate_vs_sep_s1_v;
        h->SetPoint(h->GetN(), bpm.x, bbc.rate);
        h->SetPointError(h->GetN()-1, bpm.x_err, bbc.rate_err);
        v->SetPoint(v->GetN(), bpm.y, bbc.rate);
        v->SetPointError(v->GetN()-1, bpm.y_err, bbc.rate_err);
        sum_h2 += fabs(bpm.x);
      }
    } else {
      // no bpm data
    }
  }
  if(sum_h1 > sum_h2) {
    horizontal_scan_first_ = true; 
  } else {
    horizontal_scan_first_ = false;
  }
  return 0;
}

bool BeamWidthTime::IsStep(long int time ) {
  for(auto i = steps_epoch_boundaries_.begin(); i != steps_epoch_boundaries_.end(); ++i){
    if(time >= (*i).first && time <= (*i).second) return true;  
  }
  return false;
}

int BeamWidthTime::Run() {
  CorrelateTime(); /** we know which scan is first after this is run */
  MakeScanSteps();
  FitBeamWidth();
  return 0;
}

int BeamWidthTime::MakeScanSteps(){
  int step_i = 0;
  // For each step, create a distribution of values. This distribution will be used
  // to generate an average value for the step + an RMS. Then, we can assign statistical
  // and systematic errors for each vernier scan step.
  for(auto i = steps_epoch_boundaries_.begin(); i != steps_epoch_boundaries_.end(); ++i){
    long int time_start = (*i).first;
    long int time_end   = (*i).second;
    std::vector<double> bbc_rate; // for calculation of rate fluctuation (systematic)
    std::vector<double> beam_x; // for calculation of rate fluctuation 
    std::vector<double> beam_y; // for calculation of rate fluctuation
    double bbc_gl1p = 0; // for calcualtion of rate + rate statistical error 
    double clock_gl1p = 0; // for calculation of rate + rate statistical error 

    for( long int time_i = time_start; time_i <= time_end; time_i++) {
      auto bbc_lookup = bbc_data_.find(time_i);
      if(bbc_lookup != bbc_data_.end()){  // bbc found
        bbc_rate.push_back((bbc_lookup->second).rate); // determine fluctation of rate
        bbc_gl1p   += (double)(bbc_lookup->second).bbc_gl1p;     // determine overall rate
        clock_gl1p += (double)(bbc_lookup->second).clock_gl1p; // determine overall rate
        auto bpm_lookup = bpm_data_.find(time_i);
        if( bpm_lookup != bpm_data_.end()){ //bpm found
          beam_x.push_back((bpm_lookup->second).x); // average for position, RMS for systematic
          beam_y.push_back((bpm_lookup->second).y); // average for position, RMS for systematic
        } 
      } else { // bbc not found
        std::cout << "not found bbc: " << time_i << std::endl;
      }
    }
    std::stringstream name;
    std::stringstream title;
    // Fill histograms 
    name << "bbcrate_" << step_i; 
    title << "BBC Rate Distribution for step: " << step_i << ";Rate (Hz);Counts";

    double min = *std::min_element(bbc_rate.begin(),bbc_rate.end());
    double max = *std::max_element(bbc_rate.begin(),bbc_rate.end());
    double minx = *std::min_element(beam_x.begin(), beam_x.end());
    double miny = *std::min_element(beam_y.begin(), beam_y.end());
    double maxx = *std::max_element(beam_x.begin(), beam_x.end());
    double maxy = *std::max_element(beam_y.begin(), beam_y.end());
    TH1F* bbc = new TH1F(name.str().c_str(), title.str().c_str(), 20, min-0.5*(max-min),max*2.);
    name.str(""); title.str("");
    name << "bpmx_" << step_i;
    title << "Beam Position Distribution for step: " << step_i << ";beam position (microns);Measurement";
    TH1F* bpmx = new TH1F(name.str().c_str(), title.str().c_str(), 20, minx-0.5*(maxx-minx),maxx*2.);
    name.str(""); title.str("");
    name << "bpmy_" << step_i;
    title << "Beam Position Distribution for step: " << step_i << ";beam position (microns);Measurement";
    TH1F* bpmy = new TH1F(name.str().c_str(), title.str().c_str(), 20, miny-0.5*(maxy-miny),maxy*2.);
    plot_registry_.push_back(bbc);
    plot_registry_.push_back(bpmx);
    plot_registry_.push_back(bpmy);
    for(auto rate = bbc_rate.begin(); rate!=bbc_rate.end();++rate) {
      bbc->Fill(*rate);
    }
    for(auto bpm_i = beam_x.begin(); bpm_i != beam_x.end(); ++bpm_i){
      bpmx->Fill(*bpm_i);
    }  
    for(auto bpm_i = beam_y.begin(); bpm_i != beam_y.end(); ++bpm_i){
      bpmy->Fill(*bpm_i);
    }  
    // We have now created distributions of the bbc rate, horizontal beam separation
    // and vertical beam separatio for a single vernier scan step.
    BeamWidthData bw;
    bw.rate = bbc_gl1p/clock_gl1p * CLOCK_RATE;
    bw.rate_sys_err = 0. ;/* wp - this is large for some runs and not others, ignore for now// bbc->GetRMS(); */
    double bbc_frac_err   = pow(bbc_gl1p  ,0.5)/bbc_gl1p;
    double clock_frac_err = pow(clock_gl1p,0.5)/clock_gl1p;
    bw.rate_stat_err = bw.rate*pow(pow(bbc_frac_err,2.0)+pow(clock_frac_err,2.0),0.5);
    bw.x = bpmx->GetMean();
    bw.y = bpmy->GetMean();
    bw.x_stat_err = bpmx->GetRMS();
    bw.y_stat_err = bpmy->GetRMS();
    bw.x_sys_err = 10; // from technical note
    bw.y_sys_err = 10; // from technical note

    // Now, use scan order to properly define planned scan steps 
    // NOTE that epoch_step_boundaries_ MUST be the same size as our planned steps vector.
    if( steps_epoch_boundaries_.size() != planned_beam_separation_.size() ) {
      std::cerr << " the number of planned steps does not equal the number of steps we "
        << "have extracted in the vernier scans. This is a major problem, and "
        << "this instance: " << this_name_ << " will undoubtibly crash." << std::endl;
    } else {
      if(horizontal_scan_first_) { // horizontal scan is first
        if((unsigned int)step_i < planned_beam_separation_.size()/2) { // first half
          bw.x_planned = planned_beam_separation_[step_i];
          bw.y_planned = 0.;
        } else { // second half 
          bw.y_planned = planned_beam_separation_[step_i];
          bw.x_planned = 0.;
        }
      } else {  // vertical scan is first
        if((unsigned int)step_i < planned_beam_separation_.size()/2) { // first half
          bw.y_planned = planned_beam_separation_[step_i];
          bw.x_planned = 0.;
        } else { // second half 
          bw.x_planned = planned_beam_separation_[step_i];
          bw.y_planned = 0.;
        }
      }
    }
    beam_width_data_.push_back(bw);
    step_i++;
  } // End over epoch step boundaries (dimension == number of scanning steps)
  global_output_ << "Beam Width Data" << std::endl;
  global_output_
    << "step:" << std::setw(4) << "N"
    << std::setw(11) << "bbc_rate"
    << std::setw(11) << "bbc_sys" 
    << std::setw(11) << "bbc_stat"
    << std::setw(11) << "planned_x"
    << std::setw(11) << "bpm_x" 
    << std::setw(11) << "bpm_x_sys"
    << std::setw(11) << "bpm_x_stat"
    << std::setw(11) << "planned_y"
    << std::setw(11) << "bpm_y" 
    << std::setw(11) << "bpm_y_sys"
    << std::setw(11) << "bpm_y_stat"
    << std::endl;
  step_i = 0;
  // Planned steps return to the same x_planned, y_planned coordinate, and we 
  // know this to be 0,0. ROOT does not know how to deal with this, because we set a point
  // and then set another point to the same domain value. Thus, we take the average of the
  // rates at these points and use that as the singular value, while expanding the uncertainty of
  // this point to be the combined uncertainty of both points.
  std::vector<double> planned_rate_overlap;
  std::vector<double> planned_rate_overlap_err;
  for(auto i = beam_width_data_.begin(); i != beam_width_data_.end(); ++i) {
    global_output_ << "step:" << std::setw(4) << step_i
      << std::setw(11) << (*i).rate
      << std::setw(11) << (*i).rate_sys_err 
      << std::setw(11) << (*i).rate_stat_err
      << std::setw(11) << (*i).x_planned
      << std::setw(11) << (*i).x
      << std::setw(11) << (*i).x_sys_err
      << std::setw(11) << (*i).x_stat_err
      << std::setw(11) << (*i).y_planned
      << std::setw(11) << (*i).y
      << std::setw(11) << (*i).y_sys_err
      << std::setw(11) << (*i).y_stat_err
      << std::endl;
    step_i++;
  }
  // Now, we loop over the structure to generate plots. By this time, we know which scan is first
  // and our data structure is time ordered. 
  horizontal_scan_step_ = new TGraphErrors();
  horizontal_scan_step_->SetName("horizontal_scan_step_");
  std::stringstream title;
  title << "Horizontal Scan, Run: " << run_number_ << ";separation (microns);bbc rate (Hz)";
  horizontal_scan_step_->SetTitle(title.str().c_str());

  vertical_scan_step_ = new TGraphErrors();  
  vertical_scan_step_->SetName("vertical_scan_step_");
  title.str("");
  title << "Vertical Scan, Run: " << run_number_ << ";separation (microns);bbc rate (Hz)";
  vertical_scan_step_->SetTitle(title.str().c_str());
  whole_scan_step_ = new TGraph2DErrors();     
  whole_scan_step_->SetName("whole_scan_step_");
  title.str("");
  title << "Whole Scan, Run: " << run_number_ << ";horizontal separation (microns);vertical separation (microns);bbc rate (Hz)";
  whole_scan_step_->SetTitle(title.str().c_str());

  planned_horizontal_scan_step_ = new TGraphErrors();
  planned_horizontal_scan_step_->SetName("planned_horizontal_scan_step_");
  title.str("");
  title << "Horizontal Scan, planned steps, Run: " << run_number_ << ";separation(microns);bbc rate (Hz)"; 
  planned_horizontal_scan_step_->SetTitle(title.str().c_str());

  planned_vertical_scan_step_ = new TGraphErrors();
  planned_vertical_scan_step_->SetName("planned_vertical_scan_step_");
  title.str("");
  title << "Vertical Scan, planned steps, Run: " << run_number_ << ";separation(microns);bbc rate (Hz)"; 
  planned_vertical_scan_step_->SetTitle(title.str().c_str());

  whole_scan_planned_step_ = new TGraph2DErrors();     
  whole_scan_planned_step_->SetName("whole_scan_planned_step_");
  title.str("");
  title << "Whole Scan (Planned Steps),  Run: " << run_number_ << ";planned horizontal separation (microns);planned vertical separation (microns);bbc rate (Hz)";
  whole_scan_planned_step_->SetTitle(title.str().c_str());

  plot_registry_.push_back(planned_horizontal_scan_step_);
  plot_registry_.push_back(planned_vertical_scan_step_);
  plot_registry_.push_back(horizontal_scan_step_);
  plot_registry_.push_back(vertical_scan_step_);
  plot_registry_.push_back(whole_scan_step_);
  plot_registry_.push_back(whole_scan_planned_step_);

  TGraphErrors* first_scan;
  TGraphErrors* second_scan;
  TGraphErrors* first_scan_planned;
  TGraphErrors* second_scan_planned;
  if(horizontal_scan_first_) { // horizontal scan first
    first_scan = horizontal_scan_step_;
    second_scan = vertical_scan_step_;
    first_scan_planned = planned_horizontal_scan_step_;
    second_scan_planned = planned_vertical_scan_step_;

  } else { // vertical scan first
    first_scan = vertical_scan_step_;
    second_scan = horizontal_scan_step_;
    first_scan_planned = planned_vertical_scan_step_;
    second_scan_planned = planned_horizontal_scan_step_;
  }

  for(unsigned int i = 0; i < beam_width_data_.size(); i++) {
    // split into two halves
    // first half
    BeamWidthData bw = beam_width_data_[i];
    double rate_error = pow(bw.rate_stat_err,2.0)+pow(bw.rate_sys_err,2.0);
    rate_error = pow(rate_error,0.5);
    double x_err = pow(bw.x_sys_err,2.0)+pow(bw.x_stat_err,2.0);
    x_err = pow(x_err,0.5);
    double y_err = pow(bw.y_sys_err,2.0)+pow(bw.y_stat_err,2.0);
    y_err = pow(y_err,0.5);

    double scan_1_sep = 0;
    double scan_1_sep_err = 0;
    double scan_2_sep = 0; 
    double scan_2_sep_err = 0;
    double scan_1_planned_sep = 0;
    double scan_2_planned_sep = 0;

    std::vector<double>* disp_1;
    std::vector<double>* disp_2;
    std::vector<double>* rate_1;
    std::vector<double>* rate_2;

    if(horizontal_scan_first_) {
      scan_1_sep = bw.x;
      scan_1_sep_err = x_err;
      scan_2_sep = bw.y;
      scan_2_sep_err = y_err;
      scan_1_planned_sep = bw.x_planned;
      scan_2_planned_sep = bw.y_planned;
      disp_1 = &displacement_scan_x_;
      rate_1 = &rate_scan_x_;
      disp_2 = &displacement_scan_y_;
      rate_2 = &rate_scan_y_;
    } else {
      scan_1_sep = bw.y;
      scan_1_sep_err = y_err;
      scan_2_sep = bw.x;
      scan_2_sep_err = x_err;
      scan_2_planned_sep = bw.x_planned;
      scan_1_planned_sep = bw.y_planned;
      disp_1 = &displacement_scan_y_;
      rate_1 = &rate_scan_y_;
      disp_2 = &displacement_scan_x_;
      rate_2 = &rate_scan_x_;
    }
    whole_scan_step_        ->SetPoint     (whole_scan_step_        ->GetN()   , bw.x        , bw.y        , bw.rate);
    whole_scan_step_        ->SetPointError(whole_scan_step_        ->GetN()-1 , x_err       , y_err       , rate_error);
    if( fabs(bw.x_planned - bw.y_planned) > 0.0001 ) {  // i.e. different
      whole_scan_planned_step_->SetPoint     (whole_scan_planned_step_->GetN()   , bw.x_planned, bw.y_planned, bw.rate);
      std::cout << "x,y,z: " << bw.x_planned << ", " << bw.y_planned << ", " << bw.rate << std::endl;
      whole_scan_planned_step_->SetPointError(whole_scan_planned_step_->GetN()-1 , 0           , 0           , rate_error);
    } else {  // case where we are essentially redefining the point. 
      planned_rate_overlap.push_back(bw.rate);
      planned_rate_overlap_err.push_back(bw.rate);
      // Add these points in as a single point outside the loop.
    }
    if( i < (beam_width_data_.size()/2)) {// First Scan
      first_scan->SetPoint     (first_scan->GetN()  , scan_1_sep    , bw.rate);
      first_scan->SetPointError(first_scan->GetN()-1, scan_1_sep_err, rate_error);
      first_scan_planned->SetPoint     (first_scan_planned->GetN()  , scan_1_planned_sep, bw.rate   );
      first_scan_planned->SetPointError(first_scan_planned->GetN()-1, 0.                , rate_error);
      disp_1->push_back(scan_1_sep);
      rate_1->push_back(bw.rate);
    } else { // Second Scan
      second_scan->SetPoint     (second_scan->GetN()  , scan_2_sep    , bw.rate   );
      second_scan->SetPointError(second_scan->GetN()-1, scan_2_sep_err, rate_error);
      second_scan_planned->SetPoint     (second_scan_planned->GetN()  , scan_2_planned_sep, bw.rate   );
      second_scan_planned->SetPointError(second_scan_planned->GetN()-1, 0.                , rate_error);
      disp_2->push_back(scan_2_sep);
      rate_2->push_back(bw.rate);
    }
  } // End Loop over beam_width_data_

  // Now add in problematic points to the 2D planned distribution
  // From observation of my data sets, I know the problem point is at 0,0,rate
  // this will have to be handled more carefully if someone else ever tries
  // to use this code in a different scenario...yeah right, nobody is reading this
  // or using it...
  double mean_rate = 0;
  double n = 0;
  double err = 0;
  double rms = 0;
  for(unsigned int i = 0; i < planned_rate_overlap.size(); i++) {
    mean_rate += planned_rate_overlap[i];
    n++;
    err+=pow(planned_rate_overlap_err[i],2.0); 
  }
  mean_rate/=n;
  err = pow(err,0.5);
  for(unsigned int i = 0; i < planned_rate_overlap.size(); i++){
    rms += fabs(mean_rate - planned_rate_overlap[i]);
  }
  rms/=n;
  err = pow(rms*rms + err*err,0.5);
  whole_scan_planned_step_->SetPoint(whole_scan_planned_step_->GetN(),0,0,mean_rate);
  whole_scan_planned_step_->SetPointError(whole_scan_planned_step_->GetN()-1,0,0,err);
  // wp check out if we can add the absolute steps to the above shit
  return 0;
}

int BeamWidthTime::FitBeamWidth() {
  double max_rate_x = *std::max_element(rate_scan_x_.begin(), rate_scan_x_.end());
  double min_rate_x = *std::min_element(rate_scan_x_.begin(), rate_scan_x_.end());
  double max_disp_x = *std::max_element(displacement_scan_x_.begin(), displacement_scan_x_.end()); 
  double min_disp_x = *std::min_element(displacement_scan_x_.begin(), displacement_scan_x_.end()); 
  double min_disp_plan = *std::min_element(planned_beam_separation_.begin(),planned_beam_separation_.end());
  double max_disp_plan = *std::max_element(planned_beam_separation_.begin(),planned_beam_separation_.end());

  double max_rate_y = *std::max_element(rate_scan_y_.begin(), rate_scan_y_.end());
  double min_rate_y = *std::min_element(rate_scan_y_.begin(), rate_scan_y_.end());
  double max_disp_y = *std::max_element(displacement_scan_y_.begin(), displacement_scan_y_.end());
  double min_disp_y = *std::min_element(displacement_scan_y_.begin(), displacement_scan_y_.end());

  global_output_ << "Global features of vernier scan, Run: " << run_number_ << std::endl;
  global_output_ << "  Horizontal Scan : (Min,Max) (" << min_disp_x << ", " << max_disp_x << ")" << std::endl;
  global_output_ << "  Vertical Scan   : (Min,Max) (" << min_disp_y << ", " << max_disp_y << ")" << std::endl;
  global_output_ << "  Rate (hscan)    : (Min,Max) (" << min_rate_x << ", " << max_rate_x << ")" << std::endl;
  global_output_ << "  Rate (vscan)    : (Min,Max) (" << min_rate_y << ", " << max_rate_y << ")" << std::endl;
  global_output_ << "  Planned Scan    : (Min,Max) (" << min_disp_plan << ", " << max_disp_plan << ")" << std::endl;

  fit_beam_width_x_gaus_       = new TF1("fit_beam_width_x_gaus_","gaus",min_disp_x,max_disp_x);
  fit_beam_width_y_gaus_       = new TF1("fit_beam_width_y_gaus_","gaus",min_disp_y,max_disp_y);
  fit_beam_width_plan_x_gaus_  = new TF1("fit_beam_width_plan_x_gaus_","gaus",min_disp_plan,max_disp_plan);
  fit_beam_width_plan_y_gaus_  = new TF1("fit_beam_width_plan_y_gaus_","gaus",min_disp_plan,max_disp_plan);
  std::cout << ", " << min_disp_x << ", " << max_disp_x << ", " << min_disp_y << ", " << max_disp_y << std::endl;
  fit_beam_width_xy_gaus_      = new TF2("fit_beam_width_xy_gaus_","xygaus",min_disp_x,max_disp_x,min_disp_y,max_disp_y);
  fit_beam_width_plan_xy_gaus_ = new TF2("fit_beam_width_plan_xy_gaus_","xygaus",min_disp_plan,max_disp_plan,min_disp_plan,max_disp_plan);

  double mean_x_gaus_          = (max_disp_x + min_disp_x) / 2.0;
  double mean_y_gaus_          = (max_disp_y + min_disp_y) / 2.0;
  double mean_plan_x_gaus_     = (max_disp_plan + min_disp_plan) / 2.0;
  double mean_plan_y_gaus_     = (max_disp_plan + min_disp_plan) / 2.0;
  double sigma_x_gaus_         = max_disp_x * 0.33 ;
  double sigma_y_gaus_         = max_disp_y * 0.33 ;
  double sigma_plan_x_gaus_    = max_disp_plan * 0.33;
  double sigma_plan_y_gaus_    = max_disp_plan * 0.33;
  double constant_x_gaus_      = max_rate_x;
  double constant_y_gaus_      = max_rate_y;
  double constant_plan_x_gaus_ = max_rate_x;
  double constant_plan_y_gaus_ = max_rate_y;

  fit_beam_width_x_gaus_     ->SetParameter("Mean"     , mean_x_gaus_          );
  fit_beam_width_y_gaus_     ->SetParameter("Mean"     , mean_y_gaus_          );
  fit_beam_width_plan_x_gaus_->SetParameter("Mean"     , mean_plan_x_gaus_     );
  fit_beam_width_plan_y_gaus_->SetParameter("Mean"     , mean_plan_y_gaus_     );
  fit_beam_width_x_gaus_     ->SetParameter("Sigma"    , sigma_x_gaus_         );
  fit_beam_width_y_gaus_     ->SetParameter("Sigma"    , sigma_y_gaus_         );
  fit_beam_width_plan_x_gaus_->SetParameter("Sigma"    , sigma_plan_x_gaus_    );
  fit_beam_width_plan_y_gaus_->SetParameter("Sigma"    , sigma_plan_y_gaus_    );
  fit_beam_width_x_gaus_     ->SetParameter("Constant" , constant_x_gaus_      );
  fit_beam_width_y_gaus_     ->SetParameter("Constant" , constant_y_gaus_      );
  fit_beam_width_plan_x_gaus_->SetParameter("Constant" , constant_plan_x_gaus_ );
  fit_beam_width_plan_y_gaus_->SetParameter("Constant" , constant_plan_y_gaus_ );

  fit_beam_width_x_gaus_     ->SetParLimits(1 , -1.0*sigma_x_gaus_       , sigma_x_gaus_             );
  fit_beam_width_y_gaus_     ->SetParLimits(1 , -1.0*sigma_y_gaus_       , sigma_y_gaus_             );
  fit_beam_width_plan_x_gaus_->SetParLimits(1 , -1.0*sigma_plan_x_gaus_  , sigma_plan_x_gaus_        );
  fit_beam_width_plan_y_gaus_->SetParLimits(1 , -1.0*sigma_plan_y_gaus_  , sigma_plan_y_gaus_        );
  fit_beam_width_x_gaus_     ->SetParLimits(2 , 0.1*sigma_x_gaus_        , 2.0*sigma_x_gaus_         );
  fit_beam_width_y_gaus_     ->SetParLimits(2 , 0.1*sigma_y_gaus_        , 2.0*sigma_y_gaus_         );
  fit_beam_width_plan_x_gaus_->SetParLimits(2 , 0.1*sigma_plan_x_gaus_   , 2.0*sigma_plan_x_gaus_    );
  fit_beam_width_plan_y_gaus_->SetParLimits(2 , 0.1*sigma_plan_y_gaus_   , 2.0*sigma_plan_y_gaus_    );
  fit_beam_width_x_gaus_     ->SetParLimits(0 , 0.5*constant_x_gaus_     , 1.5*constant_x_gaus_      );
  fit_beam_width_y_gaus_     ->SetParLimits(0 , 0.5*constant_y_gaus_     , 1.5*constant_y_gaus_      );
  fit_beam_width_plan_x_gaus_->SetParLimits(0 , 0.5*constant_plan_x_gaus_, 1.5*constant_plan_x_gaus_ );
  fit_beam_width_plan_y_gaus_->SetParLimits(0 , 0.5*constant_plan_y_gaus_, 1.5*constant_plan_y_gaus_ );

  // Store fit results in the shared pointer (option "S")
  std::cout << "FIT RESULTS: horizontal scan, bpm steps" << std::endl;
  TFitResultPtr x_gaus_ptr = horizontal_scan_step_->Fit(fit_beam_width_x_gaus_ ,"RS");
  std::cout << "FIT RESULTS: vertical scan, bpm steps" << std::endl;
  TFitResultPtr y_gaus_ptr = vertical_scan_step_  ->Fit(fit_beam_width_y_gaus_ ,"RS");
  std::cout << "FIT RESULTS: horizontal scan, planned steps" << std::endl;
  TFitResultPtr x_gaus_plan_ptr = planned_horizontal_scan_step_->Fit(fit_beam_width_plan_x_gaus_,"RS");
  std::cout << "FIT RESULTS: vertical scan, planned steps" << std::endl;
  TFitResultPtr y_gaus_plan_ptr = planned_vertical_scan_step_->Fit(fit_beam_width_plan_y_gaus_,"RS");

  if(x_gaus_ptr->IsValid()) {
    global_output_ << "1D gaussian fit for horizontal scan is a success!" << std::endl;
    horizontal_width_ = fit_beam_width_x_gaus_->GetParameter("Sigma");
  } else {
    global_output_ << "1D gaussian fit for horizontal scan is a failure!" << std::endl;
    std::cerr << "Cannot fit horizontal beam width slice, multidimensional fit will have bad seed values!"
      << std::endl;
    //return 1;
  }
  if(y_gaus_ptr->IsValid()) {
    global_output_ << "1D gaussian fit for vertical scan is a success!" << std::endl;
    vertical_width_ = fit_beam_width_y_gaus_->GetParameter("Sigma");
  } else {
    global_output_ << "1D gaussian fit for vertical scan is a failure!" << std::endl;
    std::cerr << "Cannot fit vertical beam width slice, multidimensional fit will have bad seed values!"
      << std::endl;
    //return 1;
  }
  if(x_gaus_plan_ptr->IsValid()) {
    global_output_ << "1D gaussian fit for horizontal scan (planned steps) is a success!" << std::endl;
    horizontal_width_ = fit_beam_width_plan_x_gaus_->GetParameter("Sigma");
  } else {
    global_output_ << "1D gaussian fit for horizontal scan (planned steps) is a failure!" << std::endl;
    std::cerr << "Cannot fit horizontal beam width slice (planned steps), multidimensional fit will have bad seed values!"
      << std::endl;
    //return 1;
  }
  if(y_gaus_plan_ptr->IsValid()) {
    global_output_ << "1D gaussian fit for vertical scan (planned steps) is a success!" << std::endl;
    vertical_width_ = fit_beam_width_plan_y_gaus_->GetParameter("Sigma");
  } else {
    global_output_ << "1D gaussian fit for vertical scan (planned steps) is a failure!" << std::endl;
    std::cerr << "Cannot fit vertical beam width slice (planned steps), multidimensional fit will have bad seed values!"
      << std::endl;
    //return 1;
  }

  // Determine inner/outer fit ranges for the beam widths 
  // Use the 1D plots as a seed for the 2D plots
  double x_mean  = fit_beam_width_x_gaus_->GetParameter("Mean");
  double x_const = fit_beam_width_x_gaus_->GetParameter("Constant");
  double x_sigma = fit_beam_width_x_gaus_->GetParameter("Sigma");
  double y_mean  = fit_beam_width_y_gaus_->GetParameter("Mean");
  double y_const = fit_beam_width_y_gaus_->GetParameter("Constant");
  double y_sigma = fit_beam_width_y_gaus_->GetParameter("Sigma");
  double avg_const = (x_const+y_const)*0.5;
  fit_beam_width_xy_gaus_->SetParameter(0, avg_const);
  fit_beam_width_xy_gaus_->SetParameter(1, x_mean   );
  fit_beam_width_xy_gaus_->SetParameter(2, x_sigma  );
  fit_beam_width_xy_gaus_->SetParameter(3, y_mean   );
  fit_beam_width_xy_gaus_->SetParameter(4, y_sigma  );
  // Using root default parameter numbers - be careful!
  fit_beam_width_xy_gaus_->SetParLimits(0, avg_const*  0.6   , avg_const * 1.6 );// 0  Constant
  fit_beam_width_xy_gaus_->SetParLimits(1, x_mean   *(-1.5)  , x_mean    * 1.5 );// 1  MeanX   
  fit_beam_width_xy_gaus_->SetParLimits(2, x_sigma  *  0.4   , x_sigma   * 1.6 );// 2  SigmaX  
  fit_beam_width_xy_gaus_->SetParLimits(3, y_mean   *(-1.5)  , y_mean    * 1.5 );// 3  MeanY   
  fit_beam_width_xy_gaus_->SetParLimits(4, y_sigma  *  0.4   , y_sigma   * 1.6 );// 4  SigmaY  
  std::cout << "FIT RESULTS: Simultaneous, Whole Scan, bpm steps" << std::endl;
  TFitResultPtr xy_gaus_ptr = whole_scan_step_    ->Fit(fit_beam_width_xy_gaus_,"RS"); 
  if(xy_gaus_ptr->IsValid()) {
    global_output_ << "2D gaussian fit is a success!" << std::endl;
    horizontal_width_ = fit_beam_width_xy_gaus_->GetParameter(2);
    vertical_width_ = fit_beam_width_xy_gaus_->GetParameter(4);
  } else {
    global_output_ << "2D gaussian fit is a failure!" << std::endl;
  }

  // Determine inner/outer fit ranges for the planned beam widths 
  // Use the 1D plots as a seed for the 2D plots
  x_mean  = fit_beam_width_plan_x_gaus_->GetParameter("Mean");
  x_const = fit_beam_width_plan_x_gaus_->GetParameter("Constant");
  x_sigma = fit_beam_width_plan_x_gaus_->GetParameter("Sigma");
  y_mean  = fit_beam_width_plan_y_gaus_->GetParameter("Mean");
  y_const = fit_beam_width_plan_y_gaus_->GetParameter("Constant");
  y_sigma = fit_beam_width_plan_y_gaus_->GetParameter("Sigma");

  avg_const = (x_const+y_const)*0.5;
  fit_beam_width_plan_xy_gaus_->SetParameter(0, avg_const);
  fit_beam_width_plan_xy_gaus_->SetParameter(1, x_mean   );
  fit_beam_width_plan_xy_gaus_->SetParameter(2, x_sigma  );
  fit_beam_width_plan_xy_gaus_->SetParameter(3, y_mean   );
  fit_beam_width_plan_xy_gaus_->SetParameter(4, y_sigma  );
  // Using root default parameter numbers - be careful!
  fit_beam_width_plan_xy_gaus_->SetParLimits(0, avg_const*  0.6   , avg_const * 1.6 );// 0  Constant
  fit_beam_width_plan_xy_gaus_->SetParLimits(1, x_mean   *(-1.5)  , x_mean    * 1.5 );// 1  MeanX   
  fit_beam_width_plan_xy_gaus_->SetParLimits(2, x_sigma  *  0.4   , x_sigma   * 1.6 );// 2  SigmaX  
  fit_beam_width_plan_xy_gaus_->SetParLimits(3, y_mean   *(-1.5)  , y_mean    * 1.5 );// 3  MeanY   
  fit_beam_width_plan_xy_gaus_->SetParLimits(4, y_sigma  *  0.4   , y_sigma   * 1.6 );// 4  SigmaY  
  std::cout << "FIT RESULTS: Simultaneous, Whole Scan, planned steps" << std::endl;
  TFitResultPtr xy_gaus_plan_ptr = whole_scan_planned_step_->Fit(fit_beam_width_plan_xy_gaus_,"RS"); 
  if(xy_gaus_ptr->IsValid()) {
    global_output_ << "2D gaussian fit is a (planned steps) success!" << std::endl;
    horizontal_width_ = fit_beam_width_plan_xy_gaus_->GetParameter(2);
    vertical_width_ = fit_beam_width_plan_xy_gaus_->GetParameter(4);
  } else {
    global_output_ << "2D gaussian fit is a (planned steps) failure!" << std::endl;
  }
  plot_registry_.push_back(fit_beam_width_x_gaus_);
  plot_registry_.push_back(fit_beam_width_y_gaus_);
  plot_registry_.push_back(fit_beam_width_xy_gaus_);
  plot_registry_.push_back(fit_beam_width_plan_x_gaus_);
  plot_registry_.push_back(fit_beam_width_plan_y_gaus_);
  plot_registry_.push_back(fit_beam_width_plan_xy_gaus_);

  // Print Global Fit Statistics (maybe to table later?)
  return 0;
}

int BeamWidthTime::SaveFigures(const std::string& figure_output_dir) {
  gErrorIgnoreLevel = kWarning;
  std::string tfile_name = figure_output_dir + "/" + run_number_ + "_BeamWidthPlots.root"; 
  std::string pdf_name   = figure_output_dir + "/" + run_number_ + "_BeamWidthPlots.pdf";
  std::string out_file_name = pdf_name + "[";
  TCanvas* booklet = new TCanvas("booklet","Scaler Rate Plots");
  TFile* root_out = new TFile(tfile_name.c_str(), "RECREATE");
  booklet->Print(out_file_name.c_str());
  for(auto plot_i = plot_registry_.begin(); plot_i != plot_registry_.end(); ++plot_i) {
    auto draw_obj = *plot_i;
    //std::string name = draw_obj->GetName();
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
  gErrorIgnoreLevel = kInfo;
  return 0;
}

int BeamWidthTime::MakePlannedStepsTable(const std::string& out_dir, const std::string& file_stub) {
  std::stringstream out_file_name;
  std::stringstream table_label  ;
  out_file_name <<  out_dir << "/" << file_stub << run_number_ << ".tex";
  table_label   <<  "tab:"  << file_stub << run_number_;

  std::ofstream out_file(out_file_name.str().c_str());

  out_file << "\\begin{table}" << std::endl;
  out_file << "\\centering" << std::endl;
  out_file << "\\begin{tabular}{c c c c c c c c}" << std::endl;
  out_file << "\\toprule" << std::endl;
  out_file << "\\textbf{$CAD_{x}$} & \\textbf{$BPM_{x}$ } & \\textbf{$\\Delta_{x}$} &\\textbf{$CAD_{y}$} & \\textbf{$BPM_{y}$} & \\textbf{$\\Delta_{y}$} & \\textbf{$BPM_{tot}$} & \\textbf{$\\Delta_{tot}$} \\\\" << std::endl;
  out_file << "($10^{-5} m$) & ($10^{-5} m$) & (\\%diff) & ($10^{-5} m$) & ($10^{-5} m$) & (\\%diff) & ($10^{-5} m$) & (\\%diff) \\\\" << std::endl;
  out_file << "\\midrule" << std::endl; 
  for(auto i = beam_width_data_.begin(); i != beam_width_data_.end(); ++i) {
    BeamWidthData d = *i; 
    float x_pct_diff;
    float y_pct_diff;
    float t_pct_diff;
    float t_sep = pow(d.x_planned*d.x_planned + d.y_planned*d.y_planned,0.5);
    std::stringstream diff_y;
    std::stringstream diff_x;
    std::stringstream diff_t;
    if(fabs(d.x_planned) > 0.01){
      x_pct_diff = fabs(fabs(d.x) - fabs(d.x_planned));
      x_pct_diff /= fabs(d.x_planned);
      x_pct_diff *= 100.;

      diff_x << " (" << std::fixed << std::setprecision(0) <<  x_pct_diff << " \\%)";
    } else { 
      diff_x << "";
    }
    if(fabs(d.y_planned) > 0.01){
      y_pct_diff = fabs(fabs(d.y) - fabs(d.y_planned));
      y_pct_diff /= fabs(d.y_planned);
      y_pct_diff *= 100.;
      diff_y << " (" << std::fixed << std::setprecision(0) << y_pct_diff << " \\%)";
    } else  {
      diff_y << "";
    }
    if(fabs(t_sep) > 0.01){
      t_pct_diff = fabs(pow(d.x*d.x+d.y*d.y,0.5) - pow(d.x_planned*d.x_planned + d.y_planned*d.y_planned,0.5));
      float denom = pow(d.x_planned*d.x_planned + d.y_planned*d.y_planned,0.5);
      t_pct_diff /= denom;
      t_pct_diff *= 100.;
      diff_t << " (" << std::fixed << std::setprecision(0) <<  t_pct_diff << " \\%)";
    } else {
      diff_t << "";
    }
    out_file 
      << std::fixed << std::setprecision(2) << d.x_planned << " & " 
      << std::fixed << std::setprecision(2) << d.x << " & " 
      << std::fixed << std::setprecision(2) << diff_x.str() << " & "
      << std::fixed << std::setprecision(2) << d.y_planned << " & " 
      << std::fixed << std::setprecision(2) << d.y << " & "
      << std::fixed << std::setprecision(2) << diff_y.str() << " & " 
      << std::fixed << std::setprecision(2) << pow(d.x*d.x + d.y*d.y,0.5) << " & " 
      << std::fixed << std::setprecision(2) << diff_t.str()
      << "\\\\" << std::endl;
  }
  out_file << "\\bottomrule" << std::endl;
  out_file << "\\end{tabular}" << std::endl;
  out_file << "\\caption{ BPM data compared to CAD data for run " << run_number_ << ". Columns are from left to right, we see the CAD planned horizontal beam displacement, the bpm-measured horizontal beam displacement, $||CAD_{x}| - |BPM_{x}||$ (to account for polarity flips in bpm), the CAD planned vertical beam displacement, the bpm-measured beam displacement, $||CAD_{y}| - |BPM_{y}||$, the total beam separation ($\\sqrt{BPM_{x}^2+BPM_{y}^2}$) and the difference between the measured total separation and the CAD planned total separation. Nominally, CAD promises to hold one beam fixed, and scan the other beam. Rows are each scan step planned and measured for the run. }" << std::endl;
  out_file << "\\label{" << table_label << "}" << std::endl;
  out_file << "\\end{table}" << std::endl;
  out_file.close();

  return 0;
}

int BeamWidthTime::SaveBeamWidthData(const std::string& out_dir) {
  std::stringstream out_file_name;
  out_file_name << out_dir << "/"  << run_number_ << "_BeamWidthData.txt";
  std::ofstream out_file(out_file_name.str().c_str());

  out_file_name.str("");
  out_file_name << out_dir << "/" << run_number_ << "_XOffsets.txt";
  std::ofstream out_xoff(out_file_name.str().c_str());

  out_file_name.str("");
  out_file_name << out_dir << "/" << run_number_ << "_YOffsets.txt";
  std::ofstream out_yoff(out_file_name.str().c_str());

  out_file_name.str("");
  out_file_name << out_dir << "/" << run_number_ << "_hWidth.txt";
  std::ofstream out_hwidth(out_file_name.str().c_str());

  out_file_name.str("");
  out_file_name << out_dir << "/" << run_number_ << "_vWidth.txt"; 
  std::ofstream out_vwidth(out_file_name.str().c_str());

  out_hwidth << "HORIZONTAL_BEAM_WIDTH " << GetHorizontalBeamWidth()/10000. << std::endl; 

  out_vwidth << "VERTICAL_BEAM_WIDTH " << GetVerticalBeamWidth()/10000. << std::endl;

  for(auto i = beam_width_data_.begin(); i != beam_width_data_.end(); ++i) {
    out_file << i->GetDataString() << std::endl;
    // desire for units of 'x-offset' to be in centimeters, but the planned
    // steps file gives us this in micrometers.
    out_xoff << "X_OFFSET "  << i->x_planned/10000. << std::endl; 
    out_yoff << "Y_OFFSET "  << i->y_planned/10000. << std::endl;
  }
  return 0;
}
