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
        ); /** Construct instance for bpm and bbc */
    int Run(); /** Creates data structures, populates graphs */
    int SaveFigures(const std::string& figure_output_dir); /** Draws, saves everything */
    int MakePlannedStepsTable(const std::string& out_dir, const std::string& file_stub);
    /** makes a table comparing planned steps to real steps, out_dir is where
     * the table is written, file_stub is the name of the file, which will be 
     * ultimately set to: file_stub << run_number_ << ".tex" */
  private:
    std::string this_name_;
    static const double CLOCK_RATE;
    std::string run_number_;
    std::vector<TObject*> plot_registry_;
    std::map < long int, BBCRateData > bbc_data_;
    std::map < long int, BeamSeparationData> bpm_data_;
    std::vector < std::pair < long int, long int> > steps_epoch_boundaries_;
    std::vector < float > planned_beam_separation_;
    std::pair<long int, long int> scan_1_;
    std::pair<long int, long int> scan_2_;

    // FITTING HELPER VARIABLES 
    std::vector<double> displacement_scan_x_; /** beam displacement for step_i during horizontal scan from data
    Used for determination of the maximum and minimum scan boundaries. */
    std::vector<double> displacement_scan_y_; /** beam displacement for step_i vertical scan from data
    Used for determination of maximum and minimum scan boundaries */
    std::vector<double> rate_scan_x_; /** rate value for step_i during horizontal scan from data, used to determine
    minimum and maximum BBC rates during the horizontal vernier scan */
    std::vector<double> rate_scan_y_; /** rate value for step_i during vertical from data, used to determine 
    minimum and maximum BBC rates durng the vertical vernier scan */
    // END FITTING HELPER VARIABLES

    std::stringstream global_output_; /** all output sent to cout is instead put here, so you get all the info after the class has run. */

    int CorrelateTime(); 
    /** Determines which scan happens first, setting horizontal_scan_first_ to true or false. 
     * Looks up BBC rate time index to find a bpm time index matching it. Generates plots which show
     * distribution of BBC rate for each time index associated with beam position. These figures are not
     * used in the analysis, because we can consider each scan step as a single measurement, where we 
     * capture statistics about the fluctuations in beam position, and rate, then generate a single
     * data point, with uncertainties for the BBC rate as a function of beam separation. */

    BeamPositionSteps bpm; /** organizes beam position data, time indexed */
    BBCRate bbc; /** organizes BBC rate data, and gl1p scaler data, time indexed */
    int FitBeamWidth(); /** Fits 1D graphs, uses these as a seed for 2D fits 
    also performs other fits using various parameterizations to obtain the beam width. */

    /** here we use these histograms to estimate the systematic error */
    std::vector<TH1F*> step_bpm_x_; /** average x separation for each beam step */
    std::vector<TH1F*> step_bpm_y_; /** average y separation for each beam step */
    std::vector<TH1F*> step_bbc_rate_; /** average rate for each step */

    int MakeScanSteps();
    /* plotting raw time correlated data has showen that there is a real fluctuation
     * in BBC rates. This will be part of the systematic error on rates. */

    bool horizontal_scan_first_; /** if true, horizontal scan was first, else, vertical scan was first */
    /** add absolute value of yvalues for rate_vs_sep plots. 
     * compares scan 0 and 1 for h and scan 0 and 1 for v
     * the largest of each comparison corresponding to when 
     * the scan is 'really' happening. */
    bool IsStep(long int time);/** checks if we're really on a step or not */

    TGraphErrors* rate_vs_sep_s0_h; /** scan 0, horizontal */
    TGraphErrors* rate_vs_sep_s0_v; /** scan 0, vertical   */
    TGraphErrors* rate_vs_sep_s1_h; /** scan 1, horizontal */
    TGraphErrors* rate_vs_sep_s1_v; /** scan 1, vertical   */

    TGraphErrors* horizontal_scan_step_;/** horizontal bbc rate vs beam displacement, plotted for each scan step */
    TGraphErrors* vertical_scan_step_;  /** vertical bbc rate vs beam displacement, plotted for each scan step */
    TGraphErrors* planned_horizontal_scan_step_; /** same as horizontal_scan_step_, but use planned beam separation */
    TGraphErrors* planned_vertical_scan_step_; /** same as vertical_scan_step_, but use planned beam separation */
    TGraph2DErrors* whole_scan_step_;   /** bbc rate vs horizontal and vertical separation, plotted for each scan step */
    TGraph2DErrors* whole_scan_planned_step_; /** same as whole_scan_step_ but used planned beam separation */
    TGraphErrors* horizontal_scan_time_; /** horizontal scan, bbc rate vs beam displacement, plotted for each time index  */
    TGraphErrors* vertical_scan_time_;   /** vertical scan, bbc rate vs beam displacement, plotted for each time index  */ 

    TH1F* step_time_; /** average time per scan step */
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
    
    std::vector<BeamWidthData> beam_width_data_;
    /** contains step-integrated beam width data */
};

#endif
