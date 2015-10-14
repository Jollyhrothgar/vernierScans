#ifndef __BBCRATESTEPS_H__
#define __BBCRATESTEPS_H__

// STL
#include <string>
#include <algorithm>
#include <vector>

// ROOT
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"

class BBCRateSteps
{
 public:
  BBCRateSteps();
  int Init (
      const std::string& RunNumber,
      const std::string& RootFileName,
      const std::string& BBCGL1PName,
      const std::string& ClockGL1PName,
      const std::string& ClockSumGL1PName,
      const std::string& StepsFileName,
      //const std::string& RateLossFileName,
      const std::string& BunchNumber = "", // leave blank as "" if bunch integrated.
      float livetimeRatio = 1.0
      );
  ~BBCRateSteps();
  int Run();
  int SetClockBBCLivetimeRatio(float livetimeRatio); // must be called before Run() if non-default livetime ratio is desired.
  int GetNumberOfSteps(); /** Returns number of steps in vernier scan */
  float GetRateStep(int i); /** Returns most accurate version of bbc rate using whole step range (Hz) */
  float GetRateStepErr(int i); /** Returns statistical error of most accurate bbc rate */ 
  float GetRateStepEventIndexCenter(int i); /** returns the center of the event-index range for a given step index */
  int GetRateStepEventIndexRange(int step_i, double& min_range, double& max_range); /** returns range of event index on step */
  int GetRateStepRateRange(int step_i, double& min_rate, double& max_rate); /** returns nominal range of bbc rate on a step */
  float GetRateStepRMS(int step_i); /* returns the rms of the BBC rate on a step */

 private:
  int Initialize();
  int LoadSteps();
  int GenerateRate();
  float livetimeRatio; // clock_livetime/bbc_livetime. Default is 1.
  int GenerateBBCRateStepHistograms(); /** BBC Rate distributions for individual steps */
  static const float clockFrequency; // Hz
  TFile* rootFile;
  TH1F* bbcGL1P;      // owned by rootFile
  TH1F* clockGL1P;    // owned by rootFile
  TH1F* clockSumGL1P; // owned by rootFile
  TH1F* bbcRate;      // owned by rootFile
  TF1* bbcRateLossCorrection;   // plug in bin number, add to BBC_RATE to adjust for losses.
  TH1F* bbcRateCorrected;
  std::vector<TH1F*> bbcRateSteps;
  std::vector<std::pair<double,double>> bbcRateRange; /** stores the range of BBC Rate values for each step */
  std::pair<int,int> prescan_range;
  std::pair<int,int> postscan_range;
  int loss_fit_start;
  int loss_fit_end;
  std::string RunNumber;
  std::string RootFileName;
  std::string BBCGL1PName;
  std::string ClockGL1PName;
  std::string ClockSumGL1PName;
  std::string StepsFileName;
  std::string RateLossFileName;
  std::string BunchNumber;
  std::vector<std::pair<int,int> > StepBoundaries;
  TCanvas* RateData;
  TCanvas* RateAdjustment;
  TCanvas* StepDistribtuions;
  std::vector<float> rateSteps;   // store bbc rate steps
  std::vector<float> rateStepsErr;// store error on bbc rate steps
  TGraphErrors* BBC_RATE; 
};
#endif
