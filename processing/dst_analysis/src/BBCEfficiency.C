#include <string>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMathBase.h"

#include "BBCEfficiency.h"
#include "VernierTreeVariables.h"

BBCEfficiency::BBCEfficiency() {}

BBCEfficiency::~BBCEfficiency() {
  // Delete everything
  for(auto i = plot_registry_.begin(); i != plot_registry_.end(); ++i) {
    if(*i) delete *i;
  }
  for(auto i = canvas_registry_.begin(); i != canvas_registry_.end(); ++i) {
    if(*i) delete *i;
  }
  if(root_file_) delete root_file_;
}

int BBCEfficiency::Init(
    const std::string& run_number_,
    const std::string& root_file_name_,
    const std::string& tree_name_
    ){
  this->run_number_     = run_number_;
  this->root_file_name_ = root_file_name_;
  this->tree_name_      = tree_name_;

  // Run 12 Defaults
  bbc_narrow_trig_bit_ = 0x00000001;
  bbc_wide_trig_bit_   = 0x00000002;
  zdc_wide_trig_bit_   = 0x00000004;

  // Default trigger range, need to update these with vertex range fit
  z_vtx_cut_min_ = -30.;
  z_vtx_cut_max_ = 30.;
  z_vtx_full_max_ = 100.;
  z_vtx_full_min_ = -100.;

  // generate map for all trigger bits to check which one fired
  // We could also use a map of names/bits populated from a databse query,
  // since we know the run number via this class' initialization. This is just used
  // to observe the difference between trigscaled, trigraw, and triglive.
  for(int i = 0; i < 32; i++){
    bit_map_[i] = ipow(2,i);
    // debug std::cout << "0x" << std::setfill('0') << std::setw(8) << std::hex << bit_map_[i] << std::endl;
  }

  float vertex_min = -300.0;
  float vertex_max = 300.0;
  int bins = 600;

  // Register Histograms
  h_bbc_wide_z_ = NULL;
  h_bbc_narrow_z_ = NULL;
  h_bbc_wide_and_narrow_z_ = NULL;
  h_zdc_wide_z_ = NULL;
  h_zdc_and_bbc_z_ = NULL;
  trig_scaled_bits_ = NULL;
  trig_raw_bits_ = NULL;
  trig_live_bits_ = NULL;
  vertex_efficiency_ = NULL;
  trigger_acceptance_ = NULL;
  trigger_acceptance_corrected_ = NULL;
  trigger_acceptance_ = NULL;
  acceptance_first_derivative_ = NULL;
  vertex_efficiency_  = NULL;
  
  std::stringstream title;
  // Initialize Histograms
  title << "BBC Wide z-vtx Distribution" << " Run: " << run_number_ << ";z-vtx (cm);counts"; 
  h_bbc_wide_z_ = new TH1F(
    "h_bbc_wide_z_",
    title.str().c_str(),
    bins,vertex_min,vertex_max);
  plot_registry_.push_back(h_bbc_wide_z_);

  title.str("");
  title << "BBC Narrow z-vtx Distribution Run: " << run_number_ << ";z-vtx (cm);counts";
  h_bbc_narrow_z_ = new TH1F(
    "h_bbc_narrow_z_",
    title.str().c_str(),
    bins,vertex_min,vertex_max);
  plot_registry_.push_back(h_bbc_narrow_z_);

  title.str("");
  title << "BBC Wide and Narrow Coincidence z-vtx Distribution" << " Run: "<< run_number_<<";z-vtx (cm);counts";
  h_bbc_wide_and_narrow_z_ = new TH1F(
    "h_bbc_wide_and_narrow_z_",
    title.str().c_str(),
    bins,vertex_min,vertex_max);
  plot_registry_.push_back(h_bbc_wide_and_narrow_z_);
  
  title.str("");
  title << "ZDC Wide z-vtx Distribution" << " Run: "<< run_number_<<";z-vtx (cm);counts";
  h_zdc_wide_z_ = new TH1F(
    "h_zdc_wide_z_",
    title.str().c_str(),
    bins,vertex_min,vertex_max);
  plot_registry_.push_back(h_zdc_wide_z_);

  title.str("");  
  title << "BBC Wide and ZDC Wide Coincidence z-vtx Distribution" << " Run: "<< run_number_<<";z-vtx (cm);counts";
  h_zdc_and_bbc_z_ = new TH1F(
    "h_zdc_and_bbc_z_",
    title.str().c_str(),
    bins,vertex_min,vertex_max);
  plot_registry_.push_back(h_zdc_and_bbc_z_);

  title.str("");
  title << "Trigger Distribution, Trig scaled, Run: " <<  run_number_ << ";Trigger Bit Number;Counts";
  trig_scaled_bits_ = new TH1F(
      "trig_scaled_bits_",
      title.str().c_str(),
      32, -0.5, 31.5);
  plot_registry_.push_back(trig_scaled_bits_);

  title.str("");
  title << "Trigger Distribution, Trig live, Run: " << run_number_ << ";Trigger Bit;Counts";
  trig_live_bits_ = new TH1F(
      "trig_live_bits_",
      title.str().c_str(),
      32, -0.5, 31.5);
  plot_registry_.push_back(trig_live_bits_);

  title.str("");
  title << "Trigger Distribution, Trig raw, Run: " << run_number_ << ";Trigger Bit;Counts";
  trig_raw_bits_ = new TH1F(
      "trig_raw_bits_",
      title.str().c_str(),
      32, -0.5, 31.5);
  plot_registry_.push_back(trig_raw_bits_);
  
  title.str("");
  title << "BBCw&ZDCw/ZDCw - Vertex Efficiency, Run: " << run_number_ << ";ZDC Z Vertex (cm); Efficiency";
  vertex_efficiency_ = new TH1F(
      "vertex_efficiency_",
      title.str().c_str(),
      bins,vertex_min,vertex_max);
  vertex_efficiency_->SetMarkerStyle(7);
  vertex_efficiency_->SetMarkerColor(kBlack);
  plot_registry_.push_back(vertex_efficiency_);

  title.str("");
  title << "BBCw&BBCnarrow/BBCw - Trigger Acceptance, Run: " << run_number_ << ";BBC Z Vertex (cm);Efficiency";
  trigger_acceptance_ = new TH1F(
      "trigger_acceptance_",
      title.str().c_str(),
      bins,vertex_min, vertex_max);
  trigger_acceptance_->SetMarkerStyle(7);
  trigger_acceptance_->SetMarkerColor(kBlack);
  plot_registry_.push_back(trigger_acceptance_);

  title.str("");
  title << "First Derivative Trigger Acceptance, Run: " << run_number_ << ";BBC Z Vertex (cm);Efficiency";
  acceptance_first_derivative_ = new TH1F(
      "acceptance_first_derivative_",
      title.str().c_str(),
      bins,vertex_min, vertex_max);
  acceptance_first_derivative_->SetMarkerStyle(7);
  acceptance_first_derivative_->SetMarkerColor(kBlack);
  plot_registry_.push_back(acceptance_first_derivative_);
  
  title.str("");
  title << "First Derivative Trigger Acceptance, Run: " << run_number_ << ";BBC Z Vertex (cm);Efficiency";
  abs_acceptance_first_derivative_ = new TH1F(
      "abs_acceptance_first_derivative_",
      title.str().c_str(),
      bins,vertex_min, vertex_max);
  abs_acceptance_first_derivative_->SetMarkerStyle(7);
  abs_acceptance_first_derivative_->SetMarkerColor(kBlack);
  plot_registry_.push_back(abs_acceptance_first_derivative_);
  
  std::cout << "Calculating BBC Efficiency Correction for Run " << run_number_ << std::endl; 

  return 0;
}

// Trigger, BBC Narrow: "BBCLL1(>0 tubes)"
int BBCEfficiency::SetBBCNarrowTrigger(int trig_bit_number ){
  bbc_narrow_trig_bit_ = trig_bit_number;
  return 0;
}
// Trigger, BBC Wide  : "BBCLL1(>0 tubes) novertex" 
int BBCEfficiency::SetBBCWideTrigger(int trig_bit_number) {
  bbc_wide_trig_bit_ = trig_bit_number;
  return 0;
}
// Trigger, ZDC Wide  : "ZDCLL1wide"
int BBCEfficiency::SetZDCWideTrigger(int trig_bit_number) {
  zdc_wide_trig_bit_ = trig_bit_number;
  return 0;
}

int BBCEfficiency::FitVertexEfficiency() {
  vertex_fit_pol2_  = new TF1("vertex_fit_pol2_","pol2", z_vtx_full_min_, z_vtx_full_max_); // Position of both bbcs
  vertex_fit_gaus_  = new TF1("vertex_fit_gaus_","gaus", z_vtx_full_min_, z_vtx_full_max_); // Position of both bbcs
  vertex_efficiency_->Fit(vertex_fit_pol2_,"R");
  vertex_efficiency_->Fit(vertex_fit_gaus_,"R+");
  plot_registry_.push_back(vertex_fit_pol2_);
  plot_registry_.push_back(vertex_fit_gaus_);
  return 0;
}

int BBCEfficiency::GetUncorrectedEfficiency() {
  /*
  int bin_min = trigger_acceptance_->GetXaxis()->FindBin(z_vtx_cut_min_);  // lower range of BBC zvtx acceptance
  int bin_max = trigger_acceptance_->GetXaxis()->FindBin(z_vtx_cut_max_);  // upper range of BBC zvtx acceptance
  double trigger_range = trigger_acceptance_->Integral(bin_min, bin_max,"width"); // bin_min should be 316
  double full_range = trigger_acceptance_->Integral("width");
  std::cout << "Uncorrected Efficiency: " << trigger_range / full_range << std::endl;
  */
  return 0;
}

/** HOW TO CALCULATE EFFICIENCY 
 * 1) Perform a fit on our acceptance plot to determine 
 *   range over which to calculate efficiency
 * 2) Integrate BBCw&BBCn plot over the range we 
 *   determined in (1) while applying the vertex 
 *   correction bin by bin. We call the result 
 *   of this integration "numerator"
 * 3) Integrate BBCw over BBC z-vertex range 
 *   applying the vertex correction bin by bin. 
 *   We call the result of this integration "denominator"
 * 4) Divide "numerator" by "denominator". 
 *   The result of this calculation is the BBC efficiency */
int BBCEfficiency::GetCorrectedEfficiency() {
  double numerator_gaus = 0;
  double numerator_pol2 = 0;
  double numerator_raw  = 0;
  double denominator_gaus = 0;
  double denominator_pol2 = 0;
  double denominator_raw  = 0;
  double uncorrected_num = 0;
  double uncorrected_denom = 0;

  for(int i = 0; i < h_bbc_wide_and_narrow_z_->GetNbinsX(); i++){ // Numerator
    double z_vtx     = h_bbc_wide_and_narrow_z_ ->GetBinCenter(i);
    double z_vtx_err = h_bbc_wide_and_narrow_z_ ->GetBinWidth(i);
    double bbc_wn_counts = h_bbc_wide_and_narrow_z_ ->GetBinContent(i);
    double bbc_w_counts  = h_bbc_wide_z_            ->GetBinContent(i);
    z_vtx_err /= 2.0; // plus/minus half of bin width -> error

    double z_correction_pol2 = vertex_fit_pol2_->Eval(z_vtx);
    double z_correction_gaus = vertex_fit_gaus_->Eval(z_vtx);
    // WARNING WARNING - we have set up in constructor that 
    // z-vertex binning is identical for bbc_wide_and_narrow_ and 
    // vertex_efficiency_, if this is changed this will break.
    double z_correction_raw = vertex_efficiency_->GetBinContent(i);
    if( fabs(z_vtx) > z_vtx_full_max_) continue;

    denominator_gaus += bbc_w_counts/z_correction_gaus;
    denominator_pol2 += bbc_w_counts/z_correction_pol2;
    denominator_raw  += bbc_w_counts/z_correction_raw;
    uncorrected_denom += bbc_w_counts;
    if( (z_vtx > z_vtx_cut_min_) && (z_vtx < z_vtx_cut_max_) ){
      numerator_gaus += bbc_wn_counts/z_correction_gaus;
      numerator_pol2 += bbc_wn_counts/z_correction_pol2;
      numerator_raw  += bbc_wn_counts/z_correction_raw;
      uncorrected_num += bbc_wn_counts;
    } 
  }

  std::cout << "BBC EFFICIENCY: " << std::endl
      << "  gaus: " << numerator_gaus  / denominator_gaus   << std::endl
      << "  pol2: " << numerator_pol2  / denominator_pol2   << std::endl
      << "  raw : " << numerator_raw   / denominator_raw    << std::endl
      << "  uncr: " << uncorrected_num / uncorrected_denom << std::endl;
  return 0;
}

int BBCEfficiency::Run(){
  // We fill the histograms by looping over the tree.
  std::cout << "Generating bbc efficiency plots...this may take a few seconds.." << std::endl;
  VernierTreeVariables v;
  TFile * f = new TFile(root_file_name_.c_str(), "READ");
  if(!f) {
    std::cout << "error opening " << root_file_name_ << " continuing to run this program will result in segfaults/errors."
        << std::endl;
    return 1;
  }
  TTree* t = (TTree*)f->Get(tree_name_.c_str());
  if(!t) {
    std::cout << "error retrieving " << tree_name_ 
      << "from " << root_file_name_ 
      << " continuing to run this program will result in segfaults/errors."
      << std::endl;
  }

  v.LinkTree(t,"READ");
  for(Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    for(auto trigger = bit_map_.begin(); trigger != bit_map_.end(); ++trigger) {
      int bit = trigger->second;
      int bit_nb = trigger->first;
      if( (bit & v.trigscaled) > 0){
        trig_scaled_bits_->Fill(bit_nb);
      }
      if( (bit & v.triglive) > 0) {
        trig_live_bits_->Fill(bit_nb);
      }
      if ((bit & v.trigraw) > 0 ) {
        trig_raw_bits_->Fill(bit_nb);
      }
      //std::cout 
      //    << "0x" << std::setfill('0') << std::setw(8) << v.trigscaled << ", "
      //    << "0x" << std::setfill('0') << std::setw(8) << v.triglive   << ", "
      //    << "0x" << std::setfill('0') << std::setw(8) << v.trigraw    << ", "
      //    << std::endl;
    }
    if( (v.trigscaled&bbc_wide_trig_bit_) > 0) {
      h_bbc_wide_z_->Fill(v.bbc_z);
      if( (v.triglive&bbc_narrow_trig_bit_) > 0 ) {
        h_bbc_wide_and_narrow_z_->Fill(v.bbc_z);
      }
    }
    if( (v.trigscaled&bbc_narrow_trig_bit_) > 0 ) {
      h_bbc_narrow_z_->Fill(v.bbc_z);
    }
    if( (v.trigscaled&zdc_wide_trig_bit_) > 0) {
      h_zdc_wide_z_->Fill(v.zdcll1_z);
      if( v.triglive&bbc_wide_trig_bit_) {
        h_zdc_and_bbc_z_->Fill(v.zdcll1_z);
      }
    }
  } // End Loop over Tree
   
  const int bins = h_bbc_wide_z_->GetNbinsX();
  // Trigger Acceptance
  for( int i = 1; i < bins+1; i++)
  {
    double bbc_wide_hits  = h_bbc_wide_z_->GetBinContent(i);
    double bbc_wide_and_narrow_hits = h_bbc_wide_and_narrow_z_->GetBinContent(i);

    if( bbc_wide_and_narrow_hits > 0 && bbc_wide_hits > 0 )
    {
      // the terms in the ratio are correlated, so to obtain error, we must take
      // covarience into account, and use the 'ratioCorrelatedError' calcualtion.
      double ratio = bbc_wide_and_narrow_hits / bbc_wide_hits ;
      double ratioCorrelatedError = pow(pow(bbc_wide_and_narrow_hits,2)/pow(bbc_wide_hits,3),0.5);
      trigger_acceptance_->SetBinContent(i, ratio);
      trigger_acceptance_->SetBinError(i, ratioCorrelatedError);
    }  
  }

  // Vertex Efficiency
  for( int i = 1; i < bins + 1; i++)
  {
    double zdc_wide_hits = h_zdc_wide_z_->GetBinContent(i); // N2
    double bbc_and_zdc_hits = h_zdc_and_bbc_z_->GetBinContent(i); // N1
    if ((zdc_wide_hits > 0.0)&&(bbc_and_zdc_hits > 0.0))
    {
      double ratio = bbc_and_zdc_hits / zdc_wide_hits;
      // correlated or uncorrelated
      //eCorr = pow(pow(N1,2)/pow(N2,3),0.5);
      double ratioUncorrelatedError = pow((bbc_and_zdc_hits/pow(zdc_wide_hits,2.0)),0.5)*pow((1.0+bbc_and_zdc_hits/zdc_wide_hits),0.5);
      // fill bbc hists
      vertex_efficiency_->SetBinContent(i,ratio); 
      vertex_efficiency_->SetBinError(i,ratioUncorrelatedError);
    }
  }
  FitVertexEfficiency();
  return 0;
}

int BBCEfficiency::FitAcceptance(float left_low = -60., float left_high = 0., float right_low = 0., float right_high = 60.) {
  // these fit ranges should encompass any vertex cut range.
  TF1* left = new TF1("left","gaus",left_low,left_high);
  left->SetLineColor(kGreen);
  TF1* right = new TF1("right","gaus",right_low,right_high);
  right->SetLineColor(kBlue);
  Double_t par[6];
  Derivative(trigger_acceptance_, acceptance_first_derivative_);

  // Needs parameter settings, etc, but we've already successfully fit.
  // left_erf_  = new TF1("left_erf_","[0]+[1]*TMath::Erf((x-[2])/[3])",-60,0);
  // right_erf_ = new TF1("right_erf_","[0]+[1]*TMath::Erf((x-[2])/[3])",0,60);

  for(int i = 0; i < acceptance_first_derivative_->GetNbinsX(); i++) {
    abs_acceptance_first_derivative_->SetBinContent(i+1,fabs(acceptance_first_derivative_->GetBinContent(i+1)));
    abs_acceptance_first_derivative_->SetBinError(i+1,acceptance_first_derivative_->GetBinError(i+1));
  }

  abs_acceptance_first_derivative_->Fit(left,"R");
  abs_acceptance_first_derivative_->Fit(right,"R+");
  left->GetParameters(&par[0]);
  right->GetParameters(&par[3]);

  // store for access later in the class to member variables.
  left_gaus_ = left;
  right_gaus_ = right;
  
  z_vtx_cut_min_ = left_gaus_->GetParameter(1); 
  z_vtx_cut_max_ = right_gaus_->GetParameter(1); 

  std::cout << "Trigger vertex cut: " << z_vtx_cut_min_ << ", " << z_vtx_cut_max_ << std::endl;

  plot_registry_.push_back(left);
  plot_registry_.push_back(right);
  return 0;
}

int BBCEfficiency::MakeFigures( std::string figure_output_dir = "/direct/phenix+spin2/beaumim/VernierAnalysis/figures/efficiencyCorrection") {
  std::string tfile_name = figure_output_dir + "/" + run_number_ + "_BBCEfficiencyPlots.root"; 
  std::string name       = figure_output_dir + "/" + run_number_ + "_BBCEfficiencyPlots.pdf";
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

int BBCEfficiency::Derivative(TH1F*& h, TH1F*& d){
  for(int i = 2; i < h->GetNbinsX()-1; i++) {
    double y1     = (double)h->GetBinContent(i-1);
    double y2     = (double)h->GetBinContent(i);
    double y1_err = (double)h->GetBinError(i-1);
    double y2_err = (double)h->GetBinError(i);
    double x1     = (double)h->GetXaxis()->GetBinCenter(i-1);
    double x2     = (double)h->GetXaxis()->GetBinCenter(i);
    // x's error is the bin width, its already taken into account
    if( fabs(x2 - x1) < 0.00001 ) {
      std::cout << "potential undefined derivative at bin: " << i << std::endl;
      std::cout << "  singularity at: " << x2 << ", " << y2 << std::endl;
      continue;
    }
    double slope = (y2-y1)/(x2-x1);
    double slope_err = pow(y1_err*y1_err + y2_err*y2_err,0.5);
    d->SetBinContent(i,slope);
    d->SetBinError(i,slope_err);
  }
  return 0;
}

int BBCEfficiency::ipow(int base, int exp) {
  if( exp == 0 ) return 1;
  if( exp == 1 ) return base;
  return base * pow(base,exp-1);
}
