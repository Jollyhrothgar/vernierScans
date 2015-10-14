#ifndef __BEAM_WIDTH_SCALERS_H__
#define __BEAM_WIDTH_SCALERS_H__

#include "PRDFBBCRate.h"
#include "/direct/phenix+spin2/beaumim/vernierScans/vernier_analysis/src/BeamPositionSteps.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TF1.h"
#include "TF2.h"

class BeamWidthScalers {
  public:
    BeamWidthScalers();
    ~BeamWidthScalers();
    int Init(
      const std::string& runNumber ,
      const std::string& bpmData   ,
      const std::string& bpmSteps  ,
      const std::string& prdf_scalers_file 
        ); /** Construct instance for bpm and bbc */
    int Run(); /** Creates data structures, populates graphs */
    int Draw(); /** Draws Graphs */
  private:
    int CorrelateTime(); /** determines which bbc time indicies belong in which step interval */
    BeamPositionSteps bpm;
    PRDFBBCRate bbc;
    int Fit(); /** Fits 1D graphs, uses these as a seed for 2D fits */
    TGraph2DErrors* beam_width_;
    TGraphErrors* first_scan_;
    TGraphErrors* second_scan_;
    TF1* fit_first_scan_;
    TF1* fit_second_scan_;
    TF2* fit_full_scan_;
    TCanvas* width_canvas_;
};

#endif
