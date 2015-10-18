#ifndef __BEAM_WIDTH_SCALERS_H__ 
#define __BEAM_WIDTH_SCALERS_H__

// DST Vernier Analysis Libraries
#include "/direct/phenix+spin2/beaumim/vernierScans/processing/dst_analysis/src/BeamPositionSteps.h"
#include "/direct/phenix+spin2/beaumim/vernierScans/processing/dst_analysis/src/BeamSeparationData.h"


// Local Analysis
#include "BBCRate.h"
#include "BBCRateData.h"
#include "BeamWidthData.h"


// ROOT Analysis
#include "TObject.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"

// STL 
#include <vector> 
#include <map> 
#include <algorithm>

class BeamWidthTime { 
 public: 
  BeamWidthTime(); 
  ~BeamWidthTime(); 
  int Init( 
      const std::string& runNumber , 
      const std::string& bpmData   ,
      const std::string& bpmSteps  , 
      const std::string& epoch_steps, 
      const std::string& prdf_scalers_file, 
      const std::string& planned_beam_steps_file
  ); 

  // After initialization is performed, filter the data into plots, load
  // various files, and do all the "work" of calculating the beam width. Run
  // is generally an organizational member, split from Init, so that any
  // conditional initialization or special options can be implemented before
  // calculating beam widths.
  int Run();

  // Loop over the plot_registry_ and save everything to the
  // figure_output_dir. Saves a root file, and pdf booklet of plots.
  int SaveFigures(const std::string& figure_output_dir);


  // makes a table comparing planned steps to real steps, out_dir is where the
  // table is written, file_stub is the name of the file, which will be
  // ultimately set to: file_stub << run_number_ << ".tex" 
  int MakePlannedStepsTable(
      const std::string& out_dir, 
      const std::string& file_stub
  );

  // Saves beam_width_data_ to a text file caled
  // "BeamWidthData_<run_number_>.txt" for use in other contexts without
  // requiring an instance of this class.
  int SaveBeamWidthData(const std::string& out_dir);

  // Returns the horizontal beam width in cm (250 microns -> 0.0250) using the
  // best model for the beam width available. Right now, this is from fitting
  // the rate vs planned separations with a xygaus.
  float GetHorizontalBeamWidth() { return horizontal_width_; };

  // Returns the vertical beam width in cm (250 microns -> 0.0250) using the
  // best model for the beam width available. Right now, this is from fitting
  // the rate vs planned separations with a xygaus.
  float GetVerticalBeamWidth() { return vertical_width_; }; 
 private: 
  std::string this_name_; 
  static const double CLOCK_RATE; 
  std::string run_number_; 
  std::vector<TObject*> plot_registry_; 
  std::map < long int, BBCRateData > bbc_data_;
  std::map < long int, BeamSeparationData> bpm_data_; 
  std::vector < std::pair < long int, long int> > steps_epoch_boundaries_; 
  std::vector < float> planned_beam_separation_; 
  std::pair<long int, long int> scan_1_;
  std::pair<long int, long int> scan_2_;

  // best possible beam width for horizontal beam (based on success/failure of
  // fits)
  float horizontal_width_;

  // best possible beam width for vertical beam (based on success/failure of
  // fits)
  float vertical_width_;

  // beam displacement for step_i during horizontal scan from data Used for
  // determination of the maximum and minimum scan boundaries. 
  std::vector<double> displacement_scan_x_; 

  // beam displacement for step_i vertical scan from data Used for
  // determination of maximum and minimum scan boundaries                 
  std::vector<double> displacement_scan_y_; 

  // rate value for step_i during horizontal scan from data, used to determine
  // minimum and maximum BBC rates during the horizontal vernier scan
  std::vector<double> rate_scan_x_; 

  // rate value for step_i during vertical from data, used to determine
  // minimum and maximum BBC rates durng the vertical vernier scan 
  std::vector<double> rate_scan_y_; 

  // This variable organizes the various diagnostic output into one streamer
  // container which can be shown at the end of execution for organizational
  // purposes.
  std::stringstream global_output_;

  // Determines which scan happens first, setting horizontal_scan_first_ to
  // true or false.  Looks up BBC rate time index to find a bpm time index
  // matching it.  Generates plots which show distribution of BBC rate for
  // each time index associated with beam position. These figures are not used
  // in the analysis, because we can consider each scan step as a single
  // measurement, where we capture statistics about the fluctuations in beam
  // position, and rate, then generate a single data point, with uncertainties
  // for the BBC rate as a function of beam separation.
  int CorrelateTime(); 

  // Organizes all beam position data, providing drawing and statistical
  // analysis of the BPM data.
  BeamPositionSteps bpm; 

  // Organizes BBC rate data, and gl1p scaler data, time indexed
  BBCRate bbc; 

  // Fits 1D graphs, uses these as a seed for 2D fits also performs other fits
  // using various parameterizations to obtain the beam width.
  int FitBeamWidth(); 

  // average x separation for each beam step
  std::vector<TH1F*> step_bpm_x_;

  // average y separation for each beam step
  std::vector<TH1F*> step_bpm_y_;

  // average rate for each step
  std::vector<TH1F*> step_bbc_rate_; 

  // plotting raw time correlated data has showen that there is a real
  // fluctuation in BBC rates. This will be part of the systematic error on
  // rates.
  int MakeScanSteps();

  // true if horizontal scan was first, false if vertical scan was first 
  bool horizontal_scan_first_;

  // searches epoch step boundaries to determine which step, if any contains
  // the epoch time "time". Returns true if one of the steps contained this
  // time. This function is used as a filter to reject points where the beam
  // may not be held fixed.
  bool IsStep(long int time);

  // first chronological scan, horizontal separation
  TGraphErrors* rate_vs_sep_s0_h;

  // first chronological scan, vertical separation
  TGraphErrors* rate_vs_sep_s0_v;

  // second chronological scan, horizontal separation
  TGraphErrors* rate_vs_sep_s1_h;

  // second chronological scan, vertical separation
  TGraphErrors* rate_vs_sep_s1_v;

  // bbc rate vs beam displacement plotted for each scan step, horizontal
  // displacements only. Beam displacement generated from bpm data
  TGraphErrors* horizontal_scan_step_;

  // bbc rate vs beam displacement plotted for each scan step, vertical
  // displacements only. Beam displacement generated from bpm data.
  TGraphErrors* vertical_scan_step_;

  // bbc rate vs planned beam displacement for each scan step, horizontal
  // displacements only. Beam displacements from CAD. 
  TGraphErrors* planned_horizontal_scan_step_;

  // bbc rate vs planned beam displacement for each scan step, vertical
  // displacements only. Beam displacements from CAD.
  TGraphErrors* planned_vertical_scan_step_; 

  // bbc rate vs horizontal and vertical displacements for all scan steps and
  // displacements, using bpm data to generate displacements.
  TGraph2DErrors* whole_scan_step_;

  // bbc rate vs horizontal and vertical displacements for all scan steos and
  // displacements, using planned displacements from CAD.
  TGraph2DErrors* whole_scan_planned_step_; 

  // bbc rate vs beam displacement during horizontal scan, plotted for each
  // time index.
  TGraphErrors* horizontal_scan_time_;

  // bbc rate vs beam displacement during vertical scan, plotted for each time
  // index.
  TGraphErrors* vertical_scan_time_;

  // Find the length of a step by taking the difference between the start of
  // the step and end of the step using the step boundaries (whch were
  // extracted from the scan step by eye). Fill this distribution with those
  // time intervals to find the average amount of time, and variation of the
  // by-eye boundaries, over a vernier scan.
  TH1F* step_time_;

  TF1* fit_beam_width_x_gaus_; 
  TF1* fit_beam_width_plan_x_gaus_; 
  TF1* fit_beam_width_x_far_gaus_; 
  TF1* fit_beam_width_x_central_gaus_;

  TF1* fit_beam_width_y_gaus_; 
  TF1* fit_beam_width_plan_y_gaus_; 
  TF1* fit_beam_width_y_far_gaus_; 
  TF1* fit_beam_width_y_central_gaus_;

  TF2* fit_beam_width_xy_gaus_; 
  TF2* fit_beam_width_plan_xy_gaus_; 
  TCanvas* width_canvas_;

  // for each scan step (an entry in the vector) we store information
  // characterizing that step. The BBC Rate and uncertainties, Beam
  // Displacements and uncertainties and planned beam displacements are all
  // included in the BeamWidthData object.
  std::vector<BeamWidthData> beam_width_data_;

};

#endif
