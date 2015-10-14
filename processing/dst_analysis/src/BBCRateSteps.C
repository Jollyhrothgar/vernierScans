// Custom
#include "BBCRateSteps.h"

// STL
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>

// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"

const float BBCRateSteps::clockFrequency = 9.5e6; // Hz
// From PRDFF I find that the rate is: 9.5281e6 -- double check!

BBCRateSteps::BBCRateSteps() {}
int BBCRateSteps::Init(
    const std::string& RunNumber,
    const std::string& RootFileName,
    const std::string& BBCGL1PName,
    const std::string& ClockGL1PName,
    const std::string& ClockSumGL1PName,
    const std::string& StepsFileName,
    //const std::string& RateLossFileName,
    const std::string& BunchNumber,
    float livetimeRatio
    )
{
  this->RunNumber = RunNumber;
  this->RootFileName = RootFileName; 
  this->BBCGL1PName = BBCGL1PName;
  this->ClockGL1PName = ClockGL1PName;
  this->ClockSumGL1PName = ClockSumGL1PName;
  this->StepsFileName = StepsFileName;
  this->RateLossFileName = RateLossFileName;
  this->BunchNumber = BunchNumber;
  this->livetimeRatio = livetimeRatio;

  std::cout << "Processing Run: " << RunNumber << ", Using root file: " << std::endl
      << "  " << RootFileName << std::endl
      << "Pulling histograms: " << std::endl
      << "  " << BBCGL1PName << std::endl
      << "  " << ClockGL1PName << std::endl
      << "  " << ClockSumGL1PName << std::endl
      << "Using Steps File: " << std::endl
      << "  " << StepsFileName << std::endl
      << "Fitting BBC Rate with range from: " << std::endl
      << "  " << RateLossFileName << std::endl
      << "  to accont for BBC Rate losses from intensity." << std::endl;
  return 0;
}

BBCRateSteps::~BBCRateSteps() {
  //delete rootFile;
  //delete BBC_RATE;
}

int BBCRateSteps::Run() {
  Initialize();
  GenerateRate();
  return 0;
}

int BBCRateSteps::GetNumberOfSteps() {
  return static_cast<int>(rateSteps.size()); 
}

float BBCRateSteps::GetRateStep(int i) {
  return rateSteps.at(i);
}

float BBCRateSteps::GetRateStepErr(int i) {
  return rateStepsErr.at(i);
}

int BBCRateSteps::SetClockBBCLivetimeRatio(float livetimeRatio) { // must be called before Run.
  if(livetimeRatio > 0) { 
    this->livetimeRatio = livetimeRatio;
  } else {
    std::cout << "Bad livetime ratio - must be greater than one. Using default ratio of 1.0" << std::endl;
    return -1;
  }
  return 0;
}

int BBCRateSteps::Initialize() {
  LoadSteps();
  rootFile = new TFile(RootFileName.c_str(),"READ");
  if(!rootFile) {
    std::cout << RootFileName << " was not located or opened correctly" << std::endl;
    return -1;
  }
  bbcGL1P      = (TH1F*) rootFile->Get(BBCGL1PName.c_str());
  if(!bbcGL1P) {
    std::cout << BBCGL1PName << " was not loaeded correctly" << std::endl;
    return -1;
  }
  clockGL1P    = (TH1F*) rootFile->Get(ClockGL1PName.c_str());
  if(!clockGL1P) {
    std::cout << ClockGL1PName << " was not loaeded correctly" << std::endl;
    return -1;
  }
  clockSumGL1P = (TH1F*) rootFile->Get(ClockSumGL1PName.c_str()); 
  if(!clockSumGL1P) {
    std::cout << ClockSumGL1PName << " was not loaeded correctly" << std::endl;
    return -1;
  }
  bbcRate = (TH1F*) rootFile->Get("BBC_RATE");
  if(!bbcRate) {
    std::cout << "BBC_RATE was not loaded correctly" << std::endl;
    return -1;
  }
  // Determine first and last filled bins for fitting purposes.
  int empty_bin = 0;
  for(int i = StepBoundaries.back().second; ; i++ ) {
    float entries = bbcRate->GetBinContent(i);
    if (fabs(entries) < 0.0001 ) {
      empty_bin = i;
      break;
    }
  }
  prescan_range.first = 1;
  prescan_range.second = StepBoundaries.front().first-100;
  postscan_range.first = StepBoundaries.back().second+100;
  postscan_range.second = empty_bin;

  RateData = new TCanvas("RateData","BBC Rate Input Data",1200,800);
  // LoadSteps was called, so we can set the draw ranges of our plots.
  bbcGL1P     ->GetXaxis()->SetRange(StepBoundaries.front().first,StepBoundaries.back().second);
  clockGL1P   ->GetXaxis()->SetRange(StepBoundaries.front().first,StepBoundaries.back().second);
  clockSumGL1P->GetXaxis()->SetRange(StepBoundaries.front().first,StepBoundaries.back().second);
  bbcRate     ->GetXaxis()->SetRange(StepBoundaries.front().first,StepBoundaries.back().second);

  RateData->Divide(3,2);
  RateData->cd(1); bbcGL1P->Draw();
  RateData->cd(2); clockGL1P->Draw();
  RateData->cd(3); clockSumGL1P->Draw();
  RateData->cd(5); bbcRate->Draw();

  std::string title = "BBC Rate, Run: " + RunNumber + " (error scaled up by factor of 75)";
  if(!BunchNumber.empty()) {
    title += " Bunch " + BunchNumber;
  }
  title += ";Vernier Scan Step # ;BBC Rate (Hz)";
  BBC_RATE = new TGraphErrors();
  BBC_RATE->SetName("BBC_RATE_GRAPH");
  BBC_RATE->SetTitle(title.c_str());
  BBC_RATE->SetMarkerStyle(7);
  BBC_RATE->SetMarkerColor(kRed);

  std::cout << "Initialization complete" << std::endl;
  return 0;
}

int BBCRateSteps::LoadSteps() {
  std::ifstream inFile(StepsFileName.c_str());
  std::string line = "";
  if(inFile) {
    while(getline(inFile,line)) {
      if(line[0] == '#')  continue;
      int step_start;
      int step_end  ;
      std::stringstream ss;
      ss.str(line);
      ss >> step_start >> step_end;
      StepBoundaries.push_back(std::make_pair(step_start,step_end));
    }
  } else {
    std::cout << "couldn't open file: " << StepsFileName << std::endl;
    return -1;
  }
  /*
  std::ifstream f(RateLossFileName.c_str());
  line = "";
  if(f) {
    while(getline(f,line)) {
      if(line[0] == '#') continue;
      loss_fit_start = 0;
      loss_fit_end = 0;
      std::stringstream ss;
      ss.str(line);
      ss >> loss_fit_start >> loss_fit_end;
    }
  } else {
    std::cout << "couldn't open file" << RateLossFileName << std::endl;
    return -1;
  }
  */
  std::cout << "Loaded " << StepBoundaries.size() << " step boundaries" << std::endl;
  //std::cout << "Will fit for loss rate over bin range: " << loss_fit_start << " to " << loss_fit_end << std::endl;
  return (int) StepBoundaries.size();
}
int BBCRateSteps::GenerateBBCRateStepHistograms() {
  std::vector<float> averages;
  std::vector<float> avg_event_index;
  TCanvas* StepDistributions = new TCanvas("StepDistributions","BBC Rate Distribution per Step",1200,800); 
  double dimension = static_cast<double>(StepBoundaries.size());
  dimension = pow(dimension,0.5) + 1.0;
  dimension = floor(dimension);
  StepDistributions->Divide(dimension,dimension);

  for(unsigned int step_i = 0; step_i < StepBoundaries.size(); step_i++) {
    int step_begin = StepBoundaries[step_i].first;
    int step_end = StepBoundaries[step_i].second;
    avg_event_index.push_back( ((double)step_end-(double)step_begin) / 2.0);
    std::stringstream name;
    std::stringstream title;

    std::vector<double> rates;
    for(int bin_i = step_begin; bin_i < step_end; bin_i++) {
      rates.push_back(bbcRate->GetBinContent(bin_i));
    }
    name << "BBCRateStep_" << step_i;
    title << "BBC Rate Step Rate Distribution, Step: " << step_i << ";BBC Rate;Count";
    auto range = std::minmax_element(rates.begin(), rates.end());
    std::pair<double,double> bbc_range = std::make_pair(*range.first, *range.second);
    bbcRateRange.push_back(bbc_range);
    TH1F* bbcRateStep = new TH1F(
        name.str().c_str(),
        title.str().c_str(), 
        10, 
        *range.first,
        *range.second
        );
    for(auto rate_entry = rates.begin(); rate_entry != rates.end(); ++rate_entry) {
      bbcRateStep->Fill(*rate_entry);
    }
    StepDistributions->cd(step_i+1);
    bbcRateStep->Draw();
    bbcRateSteps.push_back(bbcRateStep);
  }
  return 0;
}


/** This function carries out the rate loss correction
 * based on the rate loss observed from the trend in a fixed
 * range of the bbc rate as a function of event_index (event_sequence/1000)
 * We can observe rates for similar runs around each vernier scan to
 * get rate losses directly from PRDFFs, or from WCM/DCCT data.
 *
 * The rate loss is not consistent for all runs, so I have removed it.
 *
 */
int BBCRateSteps::GenerateRate() {
  int step_i = 0;

  bbcRateCorrected = (TH1F*)bbcRate->Clone("bbcRateCorrected");
  //TF1* rateLoss = new TF1("rateLoss","pol1(0)",loss_fit_start,loss_fit_end);
  //bbcRate->Fit(rateLoss,"W");
  //bbcRateLossCorrection = new TF1("bbcRateLossCorrection","pol1(0)",StepBoundaries.front().first,StepBoundaries.back().second);
  //bbcRateLossCorrection->SetParameter(0,0);
  //bbcRateLossCorrection->SetParameter(1,-1.0*rateLoss->GetParameter(1));

  bbcRate->SetLineColor(kBlue);
  //bbcRateCorrected->SetLineColor(kRed);
  //for(int i = 1; i < bbcRate->GetNbinsX(); i++) {
  //  bbcRateCorrected->SetBinContent(i,bbcRate->GetBinContent(i) + bbcRateLossCorrection->Eval((double)i));
  //}
  //RateAdjustment = new TCanvas("RateAdjustment","BBC Rate Adjustment",1200,800);
  //bbcRateCorrected->Draw();
  //bbcRate->Draw("same");
  std::cout << "========BBC RATE STEPS==========" << std::endl;
  for( auto step = StepBoundaries.begin(); step != StepBoundaries.end(); ) {
    int step_begin = step->first;
    int step_end = step->second;

    double bbc_ticks = 0;
    double clock_ticks = 0;

    // this will be ultimately more accurate than the method
    // of dealing with the BBC_RATE histogram because for certain bins
    // there will be no GL1 clock ticks (perhaps). Here, the risk is minimized
    // because we are summing over a larger time range. BBC_RATE is useful for
    // determining the boundaries of the steps.
    for(int bin = step_begin; bin <= step_end; bin++) {
      bbc_ticks += bbcGL1P->GetBinContent(bin);
      clock_ticks += clockGL1P->GetBinContent(bin);
    }
    float bbc_rate = (bbc_ticks/clock_ticks)*livetimeRatio*clockFrequency;
    float bbc_rate_err = 0;
    float bbc_err = pow(bbc_ticks,0.5);
    float clock_err = pow(clock_ticks,0.5);
    bbc_rate_err = bbc_rate*pow(pow(bbc_err/bbc_ticks,2.0)+pow(clock_err/clock_ticks,2.0),0.5);

    rateSteps.push_back(bbc_rate);
    rateStepsErr.push_back(bbc_rate_err);
    ++step;
    ++step_i;
  }
  for(unsigned int i = 0; i < rateSteps.size(); ++i) {
    std::cout << "    " << "Step " << std::setw(3) << i;
    std::cout.precision(4);
    std::cout << ", BBC Rate: " << std::setw(10) << rateSteps.at(i) << " Hz" ;
    std::cout << " +/- " << std::setw(5) << rateStepsErr.at(i) << " Hz" << std::endl;
    BBC_RATE->SetPoint(BBC_RATE->GetN(),i,rateSteps.at(i));
    float error = rateStepsErr.at(i) * 75.0;
    BBC_RATE->SetPointError(BBC_RATE->GetN()-1,0, error);
  }
  RateData->cd(4);
  BBC_RATE->Draw("AP");
  GenerateBBCRateStepHistograms();
  return 0;
}

float BBCRateSteps::GetRateStepEventIndexCenter(int i){
  float begin_step = (float) StepBoundaries[i].first;
  float end_step   = (float) StepBoundaries[i].second;
  float center = (end_step-begin_step)/2.0;
  return center;
}

int BBCRateSteps::GetRateStepEventIndexRange(int step_i, double& min_range, double& max_range){
  min_range = StepBoundaries[step_i].first;
  max_range = StepBoundaries[step_i].second;
  return 0;
}

int BBCRateSteps::GetRateStepRateRange(int step_i, double& min_rate, double& max_rate){
  min_rate = bbcRateRange[step_i].first*livetimeRatio*clockFrequency;
  max_rate = bbcRateRange[step_i].second*livetimeRatio*clockFrequency;
  return 0;
}

float BBCRateSteps::GetRateStepRMS(int step_i){
  return bbcRateSteps[step_i]->GetRMS()*livetimeRatio*clockFrequency;
}
