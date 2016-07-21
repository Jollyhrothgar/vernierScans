#include "BeamPositionSteps.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

BeamPositionSteps::BeamPositionSteps() {
  std::cout << "Creating instance of BeamPositionSteps at " << this << std::endl;
}

BeamPositionSteps::~BeamPositionSteps() {
  std::cout << "Destroying instance of BeamPositionSteps at " << this << std::endl;
  for(auto i = plot_registry_.begin(); i != plot_registry_.end(); ++i) {
    auto deleteme = *i;
    std::string name = deleteme->GetName();
    if(deleteme) {
      delete deleteme;
   //   std::cout << "Cleaning up: " << name << std::endl;
    }
  }
}

int BeamPositionSteps::Init(
    const std::string& runNumber, 
    const std::string& bpmDataFileName, 
    const std::string& stepsFileName) {
  this->runNumber = runNumber;
  this->bpmDataFileName = bpmDataFileName;
  this->stepsFileName = stepsFileName;

  std::cout << "Processing BPM Data for run: " << runNumber;
  std::cout << " Loading steps data from file:\n    " << bpmDataFileName 
      << std::endl 
      << "  and using step deliminators from:\n    " << stepsFileName 
      << std::endl;
  error_state = false;
  return 0;
}

int BeamPositionSteps::Run() {
  if (LoadBpmData()) {
    error_state = true;
    return 1; // failure returns 1
  }
  if (LoadStepsFile()) {
    error_state = true;
    return 1; // faiure returns 1
  }
  MakeSteps();
  DrawScan();
  return 0;
}

time_t BeamPositionSteps::GetScanStartTime() {
  time_t time;
  if(steps.size() == 0) return 0;
  time = steps.front().first + GetTimeOffset();
  return time;
}

time_t BeamPositionSteps::GetScanEndTime() {
  if(steps.size() == 0) return 0;
  time_t time;
  time = steps.back().second + GetTimeOffset();
  return time;
}

unsigned int BeamPositionSteps::GetNumberOfSteps() {
  return vernierScanBPMSteps.size(); 
}

int BeamPositionSteps::LoadBpmData() {
  std::cout << " Loading steps from: " << bpmDataFileName << std::endl;
  std::ifstream inFile(bpmDataFileName.c_str());
  if(inFile) {
    std::string line;
    while(getline(inFile,line)) {
      std::stringstream ss(line);
      BeamPosition bp;
      bp.Init();
      int fill = 0;
      if(line[0] == '#') continue;

      ss >> fill >> bp.time 
        >> bp.blx7 
        >> bp.bly7
        >> bp.yex7
        >> bp.yey7
        >> bp.blx8
        >> bp.bly8
        >> bp.yex8
        >> bp.yey8;
      bp.CalculateIRPositions();
      bpmData[bp.time] = bp;
      long int time = (long int)bp.time;
      beam_separation_data_[time];
      beam_separation_data_[time].Reset();
      beam_separation_data_[time].x = bp.hSeparation;
      beam_separation_data_[time].y = bp.vSeparation;
      beam_separation_data_[time].x_avg = bp.hSeparationAvg;
      beam_separation_data_[time].y_avg = bp.vSeparationAvg;
    }
  } else { 
    std::cout << "  could not open " << bpmDataFileName << std::endl;
    return 1;
  }
  return 0;
}

int BeamPositionSteps::LoadStepsFile() {
  std::ifstream inFile(stepsFileName.c_str());
  if(inFile) {
    std::string line; 
    while(getline(inFile,line)) {
      if(line[0] == '#') continue;
      std::stringstream ss(line);
      int first = 0;
      int second = 0;
      ss >> first >> second;
      steps.push_back(std::make_pair(first,second));
    }
  } else {
    std::cout << "  could not open " << stepsFileName << std::endl;
    return 1;
  }
  return 0;
}

int BeamPositionSteps::SaveEpochSteps(const std::string& out_file_directory) {
  std::string out_file_name = out_file_directory +"/" + runNumber + "_epoch_steps.txt"; 
  std::ofstream out_file(out_file_name.c_str());
  out_file << "# these steps were generated from " << stepsFileName << std::endl
      << "# and have been converted to epoch time. " << std::endl;
  for(auto i = steps.begin(); i != steps.end(); ++i) {
    long int step_begin = 0;
    long int step_end   = 0; 
    step_begin = (*i).first  + (long int)GetTimeOffset();
    step_end   = (*i).second + (long int)GetTimeOffset();
    out_file << step_begin << " " << step_end << std::endl;
  }
  out_file.close();
  return 0;
}

int BeamPositionSteps::LookupStep(double step_time) {
  int step_index = -1;
  if(steps.size() == 0 ) {
    std::cerr << "BeamPositionSteps::steps is uninitialized!!" << std::endl;
    return -1;
  } else {
    step_index = 0;
    for(auto step_itr = steps.begin(); step_itr != steps.end(); ++step_itr) {
      double step_start = double(step_itr->first) + double(GetTimeOffset());
      double step_end = double(step_itr->second) + double(GetTimeOffset());
      //std::cout << step_start << ", " << step_end << ", " << step_end - step_start << std::endl;
      if(step_time >= step_start && step_time <= step_end) {
        return step_index;
      }
      step_index++;
    }
    // If we get here, we didn't find it.
    step_index = -1;
  }
  return step_index;
}

int BeamPositionSteps::BeamPositionLookup(double time_lookup, double rejection_threshold, double& x_sep, double& y_sep) {
  // Get first difference
  double time_diff_prev = time_lookup - double(bpmData.begin()->first);
  time_t nearest_index = bpmData.begin()->first;
  for(auto bpm_i = bpmData.begin(); bpm_i != bpmData.end(); ++bpm_i) {
    double bpm_time_current = double(bpm_i->first);
    double time_diff = fabs(bpm_time_current - time_lookup);
    if( time_diff < time_diff_prev) {
      nearest_index = bpm_i->first;
      time_diff_prev = time_diff;
    }
  }

  double threshold = fabs(time_lookup - double(nearest_index));
  if ( threshold > fabs(rejection_threshold) )  {
    x_sep = -99999.;
    y_sep = -99999.;
    return 1;
  }
  
  // Now we may extract / intripolate the beam position We want to use data on
  // either side of the actual lookup value, in the case where that value is
  // not directly found.
  auto bpm_current = bpmData.find(nearest_index);
  auto bpm_first = bpm_current;
  auto bpm_second  = bpm_current;
  time_t time_lookup_value = time_t(time_lookup);
  if ( time_lookup_value == nearest_index ) {
    x_sep = bpmData[nearest_index].hSeparation;
    y_sep = bpmData[nearest_index].vSeparation;
    return 0;
  } else { 
    // Intripolate value of bpm using weighed average between two
    //bordering points.  We look at the nearest index (and therefore also, the
    //next nearest index)
    if(time_lookup_value < nearest_index ){ 
      // take nearest_index and next lower index as bounds
      bpm_first--;
    } else { // take nearest_index and next higher index as bounds
      bpm_second++;
    }

    // linear extrapolation for h and v separation
    double time_window = double(bpm_second->first) - double(bpm_first->first);
    double h_sep_second = bpm_second->second.hSeparation;
    double h_sep_first  = bpm_first ->second.hSeparation;
    double v_sep_second = bpm_second->second.vSeparation;
    double v_sep_first  = bpm_first ->second.vSeparation;
    double h_slope = (h_sep_second-h_sep_first)/(time_window);
    double v_slope = (v_sep_second-v_sep_first)/(time_window);
    double intripolated_time = time_lookup - double(bpm_first->first);
    double h_sep_intrip = h_sep_first + (intripolated_time*h_slope);
    double v_sep_intrip = v_sep_first + (intripolated_time*v_slope);
    
    x_sep = h_sep_intrip;
    y_sep = v_sep_intrip;
  }
  return 0;
}

bool BeamPositionSteps::IsFirstScan(time_t time_index){
  time_t start = GetTimeOffset() + steps[0].second;
  time_t end   = GetTimeOffset() + steps[(steps.size())/2-1].first;
  if(time_index < start || time_index > end ) return false;
  return true;
}
bool BeamPositionSteps::IsSecondScan(time_t time_index){
  time_t start = GetTimeOffset() + steps[(steps.size())/2].second;
  time_t end   = GetTimeOffset() + steps[steps.size()-1].first;
  if(time_index < start || time_index > end ) return false;
  return true;
}

time_t BeamPositionSteps::GetCentralStepTime(int step_i) {
  if(steps.size() == 0) return 0;
  std::pair<int, int> step = steps.at(step_i);
  time_t time_begin = (time_t)step.first + GetTimeOffset();
  time_t time_end = (time_t)step.second + GetTimeOffset();

  int middle = time_end - time_begin; // taking direct average overflows the time_t type.
  middle = middle / 2;
  time_t step_time = time_begin + middle;
  //std::cout << "Step Range (epoch time) step#" << step_i << ": " << time_begin << " -> " << time_end << std::endl;
  //std::cout << "Central value: " << step_time << std::endl;
  return step_time;
}

// Make any appropriate plots here as desired.
int BeamPositionSteps::DrawScan() {
  if(error_state) return 1;
  beam_separation_h = new TGraph();
  plot_registry_.push_back(beam_separation_h);
  beam_separation_h->SetName("beam_separation_h");
  beam_separation_h->SetTitle("Hoizontal Separation;time;beam separation");
  beam_separation_h->SetMarkerStyle(7);
  beam_separation_h->SetMarkerColor(kGreen+3);
  beam_separation_v   = new TGraph();
  plot_registry_.push_back(beam_separation_v);
  beam_separation_v->SetName("beam_separation_v");
  beam_separation_v->SetTitle("Horizontal Separation;time;beam separation");
  beam_separation_v->SetMarkerStyle(7);
  beam_separation_v->SetMarkerColor(kMagenta+3);
  
  time_t start_time = bpmData.begin()->first;
  std::cout << " Run: " << runNumber << std::endl 
    << " Start Time: " << start_time << std::endl;
  
  for(auto i = bpmData.begin(); i != bpmData.end(); ++i) {
    beam_separation_h->SetPoint(beam_separation_h->GetN(),i->first - start_time,i->second.hSeparation); 
    beam_separation_v->SetPoint(beam_separation_v->GetN(),i->first - start_time,i->second.vSeparation); 
  }
  TCanvas* c = new TCanvas("c","c",1200,800);
  c->Divide(2,1);
  c->cd(1);
  beam_separation_h->Draw("AP");
  c->cd(2);
  beam_separation_v->Draw("AP");

  for(auto i = steps.begin(); i != steps.end(); ++i) {
    TLine* begin = new TLine(i->first,-1000.,i->first,1000);
    begin->SetLineColor(kGreen);
    TLine* end   = new TLine(i->second,-1000.,i->second,1000);
    end->SetLineColor(kRed);
    c->cd(1);
    begin->Draw("same");
    end->Draw("same");
    c->cd(2);
    begin->Draw("same");
    end->Draw("same");
  }
  return 0;
}

int BeamPositionSteps::MakeSteps() {
  time_t start_time = bpmData.begin()->first;
  TH1F* total_linear_extrap_rms = new TH1F("total_linear_extrap_rms","Average Linear Extrapolation Rms of all BPM data",20,-3,3);
  plot_registry_.push_back(total_linear_extrap_rms);
  TH1F* average_bpm_rms = new TH1F("average_bpm_rms","Average RMS of BPM Sector 7 and 8",20,-3,3);
  plot_registry_.push_back(average_bpm_rms);

  int step_i = 0;
  for(auto i = steps.begin(); i != steps.end(); ++i) {
    // optimize this later - currently, the file is small, we can do 
    // an exhaustive search.
    std::vector<float> x_points;
    std::vector<float> y_points;
    std::vector<float> x_points_avg;
    std::vector<float> y_points_avg;
    for( auto bpm = bpmData.begin(); bpm != bpmData.end(); ++bpm) {
      time_t time = bpm->second.time - start_time;

      if( time >= i->first && time <= i->second ) {
        x_points.push_back(bpm->second.hSeparation);
        y_points.push_back(bpm->second.vSeparation);
        x_points_avg.push_back(bpm->second.hSeparationAvg);
        y_points_avg.push_back(bpm->second.vSeparationAvg);
      }
    }
    std::stringstream x_title;
    x_title << "bpm_x_step_" << step_i;
    std::stringstream y_title;
    y_title << "bpm_y_step_" << step_i;
    std::stringstream x_title_avg;
    std::stringstream y_title_avg;
    x_title_avg << x_title.str() << "_avg";
    y_title_avg << y_title.str() << "_avg";

    // Linear Extrapolation BPM Data
    TH1F* h_x = new TH1F(
        x_title.str().c_str(),
        x_title.str().c_str(),
        3,
        *std::min_element(x_points.begin(),x_points.end())-1,
        *std::max_element(x_points.begin(),x_points.end())+1
        );
    plot_registry_.push_back(h_x);
    TH1F* h_y = new TH1F(
        y_title.str().c_str(),
        y_title.str().c_str(),
        3,
        *std::min_element(y_points.begin(),y_points.end())-1,
        *std::max_element(y_points.begin(),y_points.end())+1
        );
    plot_registry_.push_back(h_y);
    
    // Average of BPM Station 7 and 8 Method
    TH1F* h_x_avg = new TH1F(
        x_title_avg.str().c_str(),
        x_title_avg.str().c_str(),
        3,
        *std::min_element(x_points_avg.begin(),x_points_avg.end())-1,
        *std::max_element(x_points_avg.begin(),x_points_avg.end())+1
        );
    plot_registry_.push_back(h_x_avg);
    TH1F* h_y_avg = new TH1F(
        y_title_avg.str().c_str(),
        y_title_avg.str().c_str(),
        3,
        *std::min_element(y_points_avg.begin(),y_points_avg.end())-1,
        *std::max_element(y_points_avg.begin(),y_points_avg.end())+1
        );
    plot_registry_.push_back(h_y_avg);

    for(unsigned int point_i = 0; point_i < x_points.size(); point_i++) {
      h_x->Fill(x_points[point_i]);
      h_y->Fill(y_points[point_i]);
      h_x_avg->Fill(x_points_avg[point_i]);
      h_y_avg->Fill(y_points_avg[point_i]);
    }

    // now we have in range steps - we can use overall average in RMS to
    // determine uncertainty in separation
    BeamStep vernierScanStep;
    vernierScanStep.xSeparation = h_x->GetMean();
    vernierScanStep.ySeparation = h_y->GetMean();
    vernierScanStep.xSeparationRMS = h_x->GetRMS();
    vernierScanStep.ySeparationRMS = h_y->GetRMS();

    // Using the average BPM measurement
    vernierScanStep.xSeparationAvg = h_x_avg->GetMean();
    vernierScanStep.ySeparationAvg = h_y_avg->GetMean();
    vernierScanStep.xSeparationAvgRMS = h_x_avg->GetRMS();
    vernierScanStep.ySeparationAvgRMS = h_y_avg->GetRMS();

    total_linear_extrap_rms->Fill(h_x->GetRMS());
    total_linear_extrap_rms->Fill(h_y->GetRMS());
    average_bpm_rms->Fill(h_x_avg->GetRMS());
    average_bpm_rms->Fill(h_y_avg->GetRMS());

    vernierScanBPMSteps.push_back(vernierScanStep);
    step_i++;
  }
  bpm_global_rms_ = total_linear_extrap_rms->GetMean();
  bpm_global_average_rms_ = average_bpm_rms->GetMean();
  total_linear_extrap_rms->Draw();
  std::cout << "BPM DATA STEPS: " << vernierScanBPMSteps.size() << std::endl;
  return 0;
}

float BeamPositionSteps::GetHStep(int i) {
  if( (i < (int) vernierScanBPMSteps.size()) && (i >= 0) ) {
    return vernierScanBPMSteps[i].GetXSeparation();
  } else {
    std::cout << "There is no step: " << i 
        << " in the range of possible steps: (0," 
        << vernierScanBPMSteps.size() << ")" << std::endl;
    return -999.9;
  }
  return 0.;
}

float BeamPositionSteps::GetVStep(int i) {
  if( (i < (int) vernierScanBPMSteps.size() )&& (i >= 0) ) {
    return vernierScanBPMSteps[i].GetYSeparation();
  } else {
    std::cout << "There is no step: " << i 
        << " in the range of possible steps: (0," 
        << vernierScanBPMSteps.size() << ")" << std::endl;
    return -999.9;
  }
  return 0.;
}

float BeamPositionSteps::GetHStepErr(int i) {
  if( (i < (int) vernierScanBPMSteps.size()) && (i >= 0) ) {
    return vernierScanBPMSteps[i].GetXSeparationError();
  } else {
    std::cout << "There is no step: " << i 
        << " in the range of possible steps: (0," 
        << vernierScanBPMSteps.size() << ")" << std::endl;
    return -999.9;
  }
  return 0.;
}

float BeamPositionSteps::GetVStepErr(int i) {
  if( (i < (int) vernierScanBPMSteps.size()) && (i >= 0) ) {
    return vernierScanBPMSteps[i].GetYSeparationError();
  } else {
    std::cout << "There is no step: " << i 
        << " in the range of possible steps: (0," 
        << vernierScanBPMSteps.size() << ")" << std::endl;
    return -999.9;
  }
  return 0.;
}

time_t BeamPositionSteps::GetTimeOffset() {
  return bpmData.begin()->first;
}

std::map<long int, BeamSeparationData>
BeamPositionSteps::GetBeamSeparationData() {
  for(auto i = beam_separation_data_.begin(); 
      i!= beam_separation_data_.end(); 
      ++i){
    BeamSeparationData& bpm = i->second;
    bpm.x_err = bpm_global_rms_;
    bpm.y_err = bpm_global_rms_;
    bpm.x_avg_err = bpm_global_average_rms_;
    bpm.y_avg_err = bpm_global_average_rms_;
  }
  return beam_separation_data_;
}

// Maybe use later for BPM data for luminosity?
int BeamPositionSteps::SaveBeamPositionData(const std::string& out_file_name){
  for(auto i = vernierScanBPMSteps.begin();
      i != vernierScanBPMSteps.end();
      ++i){
    BeamStep bs = *i;
    std::cout << bs.xSeparation 
      << ", " << bs.xSeparationRMS 
      << ", " << bs.ySeparation 
      << ", " << bs.ySeparationRMS
      << std::endl;
  }
  return 0;
}

int BeamPositionSteps::MakeFigures(const std::string& figure_output_dir) {
  std::string tfile_name = figure_output_dir + "/" + runNumber + "_BeamPositionPlots.root"; 
  TFile* root_out = new TFile(tfile_name.c_str(), "RECREATE");
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
    delete c;
  }
  int step_counter = 0;

  // Useless graphs from when I thought that maybe the difference in the
  // average BPM measurement yielded something different than the linear
  // extrapolation. It doesn't.
  TGraphErrors* g_x_linear = new TGraphErrors();
  g_x_linear->SetName("g_x_linear");
  g_x_linear->SetTitle(";step;displacement");
  TGraphErrors* g_y_linear = new TGraphErrors();
  g_y_linear->SetName("g_y_linear");
  g_y_linear->SetTitle(";step;displacement");
  TGraphErrors* g_x_avg = new TGraphErrors();
  g_x_avg->SetName("g_x_avg");
  g_x_avg->SetTitle(";step;displacement");
  TGraphErrors* g_y_avg = new TGraphErrors();
  g_y_avg->SetName("g_y_avg");
  g_y_avg->SetTitle(";step;displacement");

  for(auto i = vernierScanBPMSteps.begin(); i != vernierScanBPMSteps.end(); ++i){
    BeamStep bpm = *i;
    g_x_linear->SetPoint(
        step_counter,
        step_counter,
        bpm.xSeparation
    );
    g_y_linear->SetPoint(
        step_counter,
        step_counter,
        bpm.ySeparation
    );
    g_x_avg->SetPoint(
        step_counter,
        step_counter,
        bpm.xSeparationAvg
    );
    g_y_avg->SetPoint(
        step_counter,
        step_counter,
        bpm.ySeparationAvg
    );
    g_x_linear->SetPointError(
        step_counter,
        0,
        bpm.xSeparationRMS
    );
    g_y_linear->SetPointError(
        step_counter,
        0,
        bpm.ySeparationRMS
    );
    g_x_avg->SetPointError(
        step_counter,
        0,
        bpm.xSeparationAvgRMS
    );
    g_y_avg->SetPointError(
        step_counter,
        0,
        bpm.ySeparationAvgRMS
    );
    step_counter++;
  }
  g_x_linear->SetMarkerColor(kRed);
  g_y_linear->SetMarkerColor(kGreen);
  g_x_avg->SetMarkerColor(kRed+3);
  g_y_avg->SetMarkerColor(kGreen+3);
  g_x_linear->SetMarkerStyle(kFullCircle);
  g_y_linear->SetMarkerStyle(kFullCircle);
  g_x_avg->SetMarkerStyle(kFullCircle);
  g_y_avg->SetMarkerStyle(kFullCircle);
  g_x_linear->Write();
  g_y_linear->Write();
  g_x_avg->Write();
  g_y_avg->Write();

  root_out->Write();
  root_out->Close();
  return 0;
}
