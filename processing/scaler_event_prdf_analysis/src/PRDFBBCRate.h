#ifndef __PRDFBBCRATE_H__
#define __PRDFBBCRATE_H__

#include "TGraph.h"
#include "TObject.h"

#include <map>
#include <utility>
#include <vector>
#include <string>
#include <sstream>

/**
 *  We can  use the raw BPM data to look up which step each BBC 
 *  Rate data point belongs to!!!
 */

class PRDFBBCRate{
 public:
  PRDFBBCRate();
  ~PRDFBBCRate();
  int Init(
    const std::string& run_number, 
    const std::string& prdf_data_file
  ); /** Contains epoch_time clock bbc_w bbc_n zdc_ scalers each line */

  int Run(); /** Loop over the data structure and make the relevant plots. Data structure: prdf_data_ */
  int SaveFigures(const std::string& figure_output_dir);
 private:
  std::string run_number_;
  std::vector<TObject*> plot_registry_; 
  struct Scaler {
    // Raw Scalers
    unsigned long long bbc_w_raw_;
    unsigned long long bbc_n_raw_;
    unsigned long long zdc_w_raw_;
    unsigned long long clock_raw_;
    // Live Scalers
    unsigned long long bbc_w_live_;
    unsigned long long bbc_n_live_;
    unsigned long long zdc_w_live_;
    unsigned long long clock_live_;
    int Init() {
      bbc_w_raw_ = 0;
      bbc_n_raw_ = 0;
      zdc_w_raw_ = 0;
      clock_raw_ = 0;
      bbc_w_live_ = 0;
      bbc_n_live_ = 0;
      zdc_w_live_ = 0;
      clock_live_ = 0;
      return 0;
    }
    std::string GetScalerConsoleString() {
      std::stringstream ss;
      ss << "Clock raw: "        << clock_raw_ 
        << ", BBC Wide raw: "    << bbc_w_raw_ 
        << ", BBC Narrow raw: "  << bbc_n_raw_ 
        << ", ZDC wide raw: "    << zdc_w_raw_
        << ", Clock live: "      << clock_live_
        << ", BBC Wide live: "   << bbc_w_live_
        << ", BBC Narrow live: " << bbc_n_live_
        << ", ZDC wide live: "   << zdc_w_live_;
      return ss.str();
    }
    std::string GetScalerSerializationString() {
      std::stringstream ss;
      ss << clock_raw_ 
        << " " << bbc_w_raw_ 
        << " " << bbc_n_raw_ 
        << " " << zdc_w_raw_
        << " " << clock_live_ 
        << " " << bbc_w_live_
        << " " << bbc_n_live_ 
        << " " << zdc_w_live_;
      return ss.str();
    }
  };

  struct Livetime{
    float bbc_w;
    float bbc_n;
    float zdc_w;
    float clock;
    int Init() {
      bbc_w = 0;
      bbc_n = 0;
      zdc_w = 0;
      clock = 0;
      return 0;
    }
    std::string GetScalerConsoleString() {
      std::stringstream ss;
      ss << "Clock livetime: "        << clock 
        << ", BBC Wide livetime: "    << bbc_w 
        << ", BBC Narrow livetime: "  << bbc_n 
        << ", ZDC wide livetime: "    << zdc_w;
      return ss.str();
    }
    std::string GetScalerSerializationString() {
      std::stringstream ss;
      ss << clock 
        << " " << bbc_w 
        << " " << bbc_n 
        << " " << zdc_w;
      return ss.str();
    }
  };
  std::map<double, Livetime> livetime_;

  /** FindLivetimeBounds 
   * Give this function the target_time, the livetime coordinate on the left hand 
   * side of target_time, the livetime coordinate on the right hand side, and get
   * back out the approximate livetime. If the target_time matches a cooridnate 
   * exactly, then you get back the exact livetime.
   *
   * Returns false if NO intripolation is required because live-time exists for 
   * searched time.
   *
   * Returns true if the key was not found, and gives you the left and right hand
   * side times bordering the time of interest (epoch) so that you can 
   * intripolate with IntripolateLivetime
   */
  bool FindLivetimeBounds(double target_time, double& left_time, double& right_time );

  /** IntripolateLivetime
   * does a linear intripolation of livetime to produce best livetime value for 
   * a given pair of live-time points, and a live-time lookup index which does not
   * exist in the dataset, but is bounded on the left and right sides by real live-
   * time values. */
  float IntripolateLivetime(
      double target_time,
      double left_time,
      double left_livetime,
      double right_time,
      double right_livetime
      );

 public:
  std::map<unsigned long long, Scaler > prdf_data_;
  TGraph* bbc_raw_rate_vs_time_;
  TGraph* bbc_live_rate_vs_time_;
  std::map<std::string,TGraph*> scaler_vs_time_;

  std::vector< std::pair<double , double> > bbc_w_raw_rate_;
  std::vector< std::pair<double , double> > bbc_w_live_rate_;
  std::vector< std::pair<double , double> > bbc_n_raw_rate_;
  std::vector< std::pair<double , double> > bbc_n_live_rate_; 
};
#endif
