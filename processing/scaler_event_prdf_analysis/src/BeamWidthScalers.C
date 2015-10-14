#include "BeamWidthScalers.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include <iostream>

BeamWidthScalers::BeamWidthScalers(){
  std::cout << "Creating instance of BeamWidthScalers at " << this << std::endl;
}

BeamWidthScalers::~BeamWidthScalers() {
  std::cout << "Destroying instance of BeamwidthScalers at " << this << std::endl;
}

int BeamWidthScalers::Init(
      const std::string& runNumber ,
      const std::string& bpmData   ,
      const std::string& bpmSteps  ,
      const std::string& prdf_scalers_file 
      ){
  bbc.Init(runNumber, prdf_scalers_file);
  bbc.Run();
  bpm.Init(
    runNumber,
    bpmData,
    bpmSteps
    );
  bpm.Run();
  
  // Maybe fastest way is to lookup?
  width_canvas_ = new TCanvas("width_canvas_", "Beam Width", 1600, 600);
  width_canvas_->Divide(3,1);
  beam_width_ = new TGraph2DErrors();
  beam_width_->SetName("beam_width_");
  beam_width_->SetTitle("Beam Width;Horizontal Separation;Vertical Separation;BBC Rate");
  beam_width_->SetMarkerStyle(kFullCircle);
  beam_width_->SetMarkerSize(0.5);
  // std::string xygaus_formula = "[0]*exp(-0.5*((x-[1])/[2])**2)*exp(-0.5*((y-[3])/[4])**2)";
  fit_full_scan_ = new TF2("fit_full_scan_","xygaus",-800,800,-800,800);

  first_scan_ = new TGraphErrors();
  second_scan_ = new TGraphErrors();

  first_scan_->SetName("first_scan_");
  second_scan_->SetName("second_scan_");

  first_scan_->SetTitle("Horizontal scan;Horizontal Separation;BBC Rate");
  second_scan_->SetTitle("Vertical scan;Vertical Separation;BBC Rate");

  first_scan_->SetMarkerStyle(kFullCircle);
  second_scan_->SetMarkerStyle(kFullCircle);

  first_scan_->SetMarkerColor(kBlue);
  second_scan_->SetMarkerColor(kGreen-2);

  fit_first_scan_ = new TF1("fit_first_scan_","gaus",-300,300);
  fit_second_scan_ = new TF1("fit_second_scan_","gaus",-300,300);

  return 0;
}

int BeamWidthScalers::CorrelateTime() {
  for( auto time_itr = bbc.bbc_n_live_rate_.begin(); time_itr != bbc.bbc_n_live_rate_.end(); ++time_itr) {
    time_t time_index = time_itr->first;
    double time_index_d = double(time_index);
    //std::cout << time_index << ", at step: " << bpm.LookupStep(time_index_d) << std::endl;
    double x_sep, y_sep;
    if( bpm.BeamPositionLookup(time_index_d,3,x_sep,y_sep) == 0 ) {// success
      double bbc_rate = time_itr->second;
      if( bpm.IsFirstScan(time_index) || bpm.IsSecondScan(time_index)) {
        beam_width_->SetPoint(beam_width_->GetN(), x_sep, y_sep, bbc_rate);
        beam_width_->SetPointError(beam_width_->GetN()-1, 10, 10, bbc_rate*0.03);
      }
      if(bpm.IsFirstScan(time_index)) {
        first_scan_->SetPoint(first_scan_->GetN(), x_sep, bbc_rate);
        first_scan_->SetPointError(first_scan_->GetN()-1, 10, bbc_rate*0.03); // BPM uncertainty is 10 microns, using 3% error for BBC 
      }
      if(bpm.IsSecondScan(time_index)){
        second_scan_->SetPoint(second_scan_->GetN(), y_sep, bbc_rate);
        second_scan_->SetPointError(second_scan_->GetN()-1, 10, bbc_rate*0.03);
      }
    }
  }
  return 0;
}

int BeamWidthScalers::Fit() {
  first_scan_->Fit(fit_first_scan_);
  second_scan_->Fit(fit_second_scan_);
  // ROOT Default Parameter Numbers for gaus
  // 0 Constant
  // 1 Mean
  // 2 Sigma
  
  // [0]*exp(-0.5*((x-[1])/[2])**2)*exp(-0.5*((y-[3])/[4])**2) 
  double first_sigma    = fit_first_scan_ ->GetParameter( 2 );
  double second_sigma   = fit_second_scan_->GetParameter( 2 );
  double first_mean     = fit_first_scan_ ->GetParameter( 1 );
  double second_mean    = fit_second_scan_->GetParameter( 1 );
  double first_constant = fit_first_scan_ ->GetParameter( 0 );
  std::cout << " first_sigma    " << first_sigma    << std::endl;
  std::cout << " second_sigma   " << second_sigma   << std::endl;
  std::cout << " first_mean     " << first_mean     << std::endl;
  std::cout << " second_mean    " << second_mean    << std::endl;
  std::cout << " first_constant " << first_constant << std::endl;
  
  // 0  Constant
  // 1  MeanX   
  // 2  SigmaX  
  // 3  MeanY   
  // 4  SigmaY  
  double mean_tolerance  = 100.0;
  double sigma_tolerance = 100.0;
  double const_tolerance = 100.0;
  fit_full_scan_->SetParLimits( 0 ,  first_constant - first_constant*const_tolerance , first_constant + first_constant*const_tolerance );
  fit_full_scan_->SetParLimits( 1 ,  first_mean - mean_tolerance*first_mean          , first_mean + mean_tolerance*first_mean );
  fit_full_scan_->SetParLimits( 2 ,  first_sigma - first_sigma*sigma_tolerance       , first_sigma + first_sigma*sigma_tolerance );
  fit_full_scan_->SetParLimits( 3 ,  second_mean - mean_tolerance*second_mean        , second_mean + mean_tolerance*second_mean);
  fit_full_scan_->SetParLimits( 4 ,  second_sigma - second_sigma*sigma_tolerance     , second_sigma + second_sigma*sigma_tolerance);
  beam_width_->Fit(fit_full_scan_,"M");
  TFile* f = new TFile("test.root","RECREATE");
  f->cd();
  fit_full_scan_->Write();
  beam_width_->Write();
  f->Write();
  f->Close();
  delete f;


  return 0;
}

int BeamWidthScalers::Draw() {
  width_canvas_->cd(1);
  gStyle->SetOptFit();
  beam_width_->Draw("APcolz");
  fit_full_scan_->Draw("same surf1z");
  width_canvas_->cd(2);
  gStyle->SetOptFit();
  first_scan_->Draw("AP");
  width_canvas_->cd(3);
  gStyle->SetOptFit();
  second_scan_->Draw("AP");
  return 0;
}

int BeamWidthScalers::Run() {
  CorrelateTime();
  Fit();
  return 0;
}
