#ifndef __BEAMWIDTH_H__
#define __BEAMWIDTH_H__

// STL
#include <string>

// ROOT
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"

// CUSTOM
#include "BBCRateSteps.h"
#include "BeamPositionSteps.h"

class BeamWidth
{
 public:
  BeamWidth();
  ~BeamWidth();
  int Init(
      const std::string& runNumber       ,
      const std::string& RootFileName    ,
      const std::string& BBCGL1PName     ,
      const std::string& ClockSumGL1PName,
      const std::string& ClockGL1PName   ,
      const std::string& StepsFileName   ,
      //const std::string& RateLossFileName,
      const std::string& BunchNumber     ,
      float livetimeRatio                ,
      const std::string& bpmData         ,
      const std::string& bpmSteps   
      );
  int Run();
 private:
  BeamPositionSteps bpm;
  BBCRateSteps bbc;
  TGraphErrors *HScan1; // Scan1: First half of vernier scan, Scan2: second half.
  TGraphErrors *HScan2;
  TGraphErrors *VScan1;
  TGraphErrors *VScan2;
  TGraphErrors *hScan; // The 'actual' horizontal scan (points to either HScan1 or HScan2)
  TGraphErrors *vScan; // The 'actual' vertical scan (points to eihter VScan1 or VScan2)
  double hScanMean     ;
  double hScanConstant ;
  double hScanSigma    ;
  double vScanMean     ;
  double vScanConstant ;
  double vScanSigma    ;
  double ScanConstant ;
  double ScanSigmaX   ;
  double ScanSigmaY   ;
  bool HorizontalScanFirst;
  TGraph2DErrors * BeamWidth2D;
  TF1* hScanFit1D;
  TF1* vScanFit1D;
  TF2* ScanFit2D;
  int MakeBeamWidth();
  int MakeCorrectedBeamWidth();
  int CorrectForRateLosses();
  int DrawBeamWidth();
  int NumberOfSteps;
  TCanvas* BeamWidthCanvas1D;
  TCanvas* BeamWidthCanvas2D;
  TCanvas* BeamWidthAll;
  // Local Data for beam correction:
  std::vector<float> rate;  
  std::vector<float> rate_step_stat_err;
  std::vector<double> rate_max; 
  std::vector<double> rate_min;
  std::vector<double> rate_ievent_min_range;
  std::vector<double> rate_ievent_max_range;
  std::vector<double> rate_ievent_center;
  std::vector<float> rate_step_rms;
  std::vector<float> hsep;
  std::vector<float> vsep;
  std::vector<float> hsep_err;
  std::vector<float> vsep_err;
};
#endif
