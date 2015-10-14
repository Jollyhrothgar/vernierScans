#ifndef __BBCRATE_H__
#define __BBCRATE_H__

#include "TGraph.h"
#include "TH1.h"
#include "TH1F.h"

#include <map>
#include <utility>
#include <vector>
#include <string>
#include <sstream>

#include "BBCRateData.h"

/**
 *  We can  use the raw BPM data to look up which step each BBC 
 *  Rate data point belongs to!!!
 */


class BBCRate{
 public:
  BBCRate();
  ~BBCRate();
  int Init(const std::string& gl1p_scalers_data_file, const std::string& run_number); 
  /** Contains epoch_time clock bbc_w bbc_n zdc_ scalers each line */
  
  int Run(); 
  /** Loop over the data structure and make the relevant plots. 
   * Data structure: prdf_data_ */

  int MakeFigures(const std::string& out_file_stub); 
  /** Draw the plots that we make in Run. Sends to root 
   * file and pdf by appending whatever you give out_file_stub 
   * with ".root" and ".pdf" */

 private:
  static const double CLOCK_RATE;
  double bbc_max_rate_;
  double bbc_min_rate_;
  long int time_offset_; /* Digit to offset time so that ROOT may plot without barfing */
  int MakeHistograms(); /* dump contents of gl1p_scalers_ to relevant histograms */

  double GetGL1PRate(long int gl1p_scaler, long int gl1p_clock); 
  /** gl1p scalers are counts, so we take the ratio 
   * of scaler to clock, multiplying by clock frequency */

  double GetGL1PRateError(long int gl1p_scaler, long int gl1p_clock);
  /** gl1p scalers are counts, so we propagate errors for  
   * counters in the standard way. */

  std::string run_number_;
  struct GL1P {      // Run 12 Settings
    long int bbc_;   // "BBCLL1(>0 tubes)"
    long int zdc_n_; // "ZDCLL1wide"              
    long int zdc_w_; // "ZDCLL1narrow"     
    long int clock_; // "Clock"     

    int Init() {
      bbc_   = 0;
      zdc_n_ = 0;
      zdc_w_ = 0;
      clock_ = 0;
      return 0;
    }
    std::string GetScalerConsoleString() {
      std::stringstream ss;
      ss << "Clock : "            << clock_ 
        << ", BBCL1(>0 tubes): "  << bbc_ 
        << ", ZDC Narrow : "      << zdc_n_
        << ", ZDC wide : "        << zdc_w_;
      return ss.str();
    }
    std::string GetScalerSerializationString() {
      std::stringstream ss;
      ss << clock_ 
        << " " << bbc_ 
        << " " << zdc_n_ 
        << " " << zdc_w_;
      return ss.str();
    }
  };
 public:
  std::map<unsigned long, GL1P > gl1p_scalers_;
  std::map<std::string, TH1*> hist;
  TH1F* bbc_rate_;
  TH1F* zdc_w_rate_;
  TH1F* zdc_n_rate_;
  TH1F* clock_rate_;
  std::vector<TObject*> plot_registry_;
  std::map<long int , BBCRateData > bbc_rate_data_;
  std::map<long int, BBCRateData> GetBBCRateData() { return bbc_rate_data_; };
};
#endif
