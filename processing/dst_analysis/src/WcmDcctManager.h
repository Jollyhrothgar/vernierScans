#ifndef __WCMDCCTMANAGER_H__
#define __WCMDCCTMANAGER_H__

// STL
#include <map>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>
#include <cstdio>
#include <fstream>
#include <sstream>

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TTree.h"
#include "TF1.h"


/*
 * This class assumes that the time stamps in any of the files used have
 * already been converted to the correct epoch time. 
 */
class WcmDcctManager{
 public:
  WcmDcctManager();
  int PrintDataTimeIndex(time_t time);
  void PrintBlueDCCT();
  void PrintYellowDCCT();
  void PrintBlueWcmTotal();
  void PrintYellowWcmTotal();
  bool CheckIfTimeExists(time_t time); 
  int Init( 
    const std::string& run_number_, 
    const std::string& wcm_b_fname_, 
    const std::string& wcm_y_fname_, 
    const std::string& dcct_fname_
  );
  int Run(); /** Generates plots from loaded data. Gets data ready to use in vernier luminosity calculation */
  int SaveFigures(const std::string& figure_output_dir); /** saves a root file and pdf of objects in rsave_registry_ to figure_output_dir */
  int SetTimeOffset(time_t offset); /* add this number to every time in the dataset, everywhere 
                                     * external option so that WCM/DCCT data may be synched up to
                                     * the bpm data. */
  int SetScanRange(time_t begin, time_t end); /** sets the epoch time range over which the vernier scan takes place */
  double GetLuminosityLoss(); /** returns luminosity_loss_ */
  int ShowSummary(); /** prints out summary for average bunch information. */
  ~WcmDcctManager();
 private:
  std::string this_name_; /** keep track of a specific instance of this class */
  std::vector<TObject*> save_registry_;
  std::vector<TObject*> bunch_registry_;

  int GenerateBounds(); /** loop over data for general characteristics of run. called in ::Run, so that you can
                          initialize with a time boundary which differs from the default */
  int LoadDcct();   
  int LoadWcmBlue();
  int LoadWcmYellow();  
  int InitPlots(const int beam_choice); /** Called in ::Run so that optional initialization can take place after
                                          required initialization */
  int MakePlots(const int beam_choice); /** Called in ::Run so that optional initialization can take place after
                                          required initialization */
  std::string run_number_;
  std::string wcm_b_fname_;
  std::string wcm_y_fname_;
  std::string dcct_fname_;
  time_t BEGIN_SCAN; /** epoch time stamp corresponding to the beginning of the vernier scan */
  time_t END_SCAN; /** epoch time stamp corresponding to the end of the vernier scan */

  int CalculateLuminosityLoss();
  /** returns 0 if successful. Uses rate_loss_ fit, BEGIN_SCAN and END_SCAN
   * to calculate how much luminosity is lost due to beam losses over the course
   * of the vernier scan. This will go in to the overall systematic error of the luminsoity. */
  double luminosity_loss_;
  static const int NUMBER_OF_BUNCHES;

  struct BeamPopulation {
    BeamPopulation() {
      std::stringstream ss; 
      ss << "BeamPopulation_0x" << std::hex << this;
      this_name_ = ss.str();
    }
    ~BeamPopulation(){};
    std::vector<float> wcm_; /** population of each bunch in beam (x10^9) */
    std::vector<float> wcm_calib_; /** calibrated population of each bunch (x10^9) */
    float dcct_; /** total beam current (x10^11) */
    float calib_; /** wcm calibration constant */
    float wcm_tot_; /** sum of wcm_ (x10^11) */
    float wcm_tot_calib_; /* calibrated sum of wcm_ (x10^11) */
    std::string this_name_;
    bool WcmAvailable() {
      if( wcm_.size() == 120) {
        return true;
      }
      return false;
    }
    float GetWcmBunch(int bunch_number) {
      return wcm_.at(bunch_number);
    }
    float GetCalibWcmBunch(int bunch_number) {
      return wcm_.at(bunch_number) * calib_;
    }
    int Reset () {
      wcm_.clear();
      dcct_ = -1;
      calib_ = -1;
      wcm_tot_ = -1;
      wcm_tot_calib_ = -1;
      return 0;
    }
    int SumWcm() {
      float sum = 0.;
      if(WcmAvailable()) {
        for(auto i = wcm_.begin(); i != wcm_.end(); ++i) {
          sum += *i;
        }
        wcm_tot_ = sum/100.0; // For promised order of magnetude
        return 0;  
      }
      return 1;
    }            
    int Calibrate(){
      calib_ = dcct_/wcm_tot_; // multiply wcm_ by this constant
      float sum = 0;
      for(int i = 0; i < 120; i++) {
        wcm_calib_.push_back(wcm_.at(i)*calib_);
	sum += wcm_calib_.back();
      } 
      wcm_tot_calib_ = sum/100.0; // For promised order of magnetude
      return 0;
    }
  };

#ifndef __CINT__
  typedef std::map<time_t, BeamPopulation> WcmDcctData;
#else
  class WcmDcctData;
#endif 

  /** data storage for the WCM and DCCT data */
  std::map<time_t, BeamPopulation> blue_pop_; /** population data blue beam */
  std::map<time_t, BeamPopulation> yell_pop_; /** population data yellow beam */

  std::map<int,TH1F*>  wcm_dist_blue_; /** calibrated wcm population for each bunch */
  std::map<int,TH1F*>  wcm_dist_yell_; /** calibrated wcm population for each bunch */
  std::map<int,TGraph*> wcm_dist_vs_time_blue_; /** each bunch population (calibrated) plotted vs time */
  std::map<int,TGraph*> wcm_dist_vs_time_yell_; /** each bunch population (calibrated) plotted vs time */ 
  TGraph* calib_vs_time_blue_; /** calibration constant vs time */
  TGraph* calib_vs_time_yell_; /** calibration constant vs time */
  TH1F* calib_blue_; /** distribution of calibration constant */
  TH1F* calib_yell_; /** distribution of calibration constant */
  TGraph* dcct_vs_time_blue_; /** DCCT population as a function of time */
  TGraph* dcct_vs_time_yell_; /** DCCT population as a function of time */
  TGraph* wcm_tot_vs_time_blue_; /** Summed WCM vs time, calibrated */
  TGraph* wcm_tot_vs_time_yell_; /** Summed WCM vs time, calibrated */
  TGraph* wcm_product_vs_time_; /* Product of Blue WCM and Yellow WCM vs time, calibrated */
  TF1* rate_loss_; /* Fit to wcm_product_vs_time_, used to determine BBC rate losses */
  int FitRateLoss(); /* Performs the fit rate_loss_ to wcm_product_vs_time_ */

  /** global features of wcm/dcct data set */
  float wcm_blue_max_; /** the maximum calibrated wcm bunch population of any single bunch in blue beam */
  float wcm_yell_max_; /** the maximum calibrated wcm bunch population of any single bunch in yellow beam */
  float dcct_blue_max_; /** the maximum dcct population blue beam */
  float dcct_yell_max_; /** the maximum dcct population yellow beam */
  float wcm_blue_min_; /** the minimum calibrated wcm bunch population of any single bunch in blue beam */
  float wcm_yell_min_; /** the minimum calibrated wcm bunch population of any single bunch in yellow beam */
  float dcct_blue_min_; /** the minimum dcct population blue beam */
  float dcct_yell_min_; /** the minimum dcct population yellow beam */
  float calib_blue_min_; /** average calibration constant of blue beam */
  float calib_blue_max_; /** rms of average calibration constant of blue beam */
  float calib_yell_min_; /** average calibration constant of yell beam */
  float calib_yell_max_; /** rms of average calibration constant of yell beam */

  time_t TIME_OFFSET;
};

#endif
