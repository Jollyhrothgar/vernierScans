#ifndef __BBC_EFFICIENCY_H__
#define __BBC_EFFICIENCY_H__

#include <string>
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <map>

#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"

class BBCEfficiency
{
 public:
  BBCEfficiency();
  ~BBCEfficiency();
  int Init(
    const std::string& run_number_,
    const std::string& root_file_name_,
    const std::string& tree_name_);

  int SetBBCNarrowTrigger(int trig_bit_number );
  int SetBBCWideTrigger(int trig_bit_number)   ;
  int SetZDCWideTrigger(int trig_bit_number)   ;
  /** Override the default trigger bits in the ctor */
  int Run();
  /** Loop over the tree, fill histograms, then call the fit */

  int MakeFigures(std::string figure_output_dir);
  /** Save all figures, create a registry of all figures generated so they can be loaded again. */

  int FitAcceptance(float left_low, float left_hi, float right_low, float right_high);
  /** Because the acceptance takes the form of a two sided error function, 
   * and because the error function depends on the integral of a gaussian
   * we can directly access this via the derivative of the trigger acceptance.
   * we take the derivative of the trigger acceptance (after smoothing) and 
   * choosing binning so that fluctions are minimized. We then fit both gaussian
   * peaks in the derivative, individually, and then simultaneously, and can use the
   * mean of the two gaussians as trigger cuts. 
   *
   * left_low <-> left_high defines range for left acceptance fit range
   * right_low <-> right_high defines range for right acceptance fit range */

  int GetCorrectedEfficiency();
  /** Gets the bbc efficiency without the zvertex correction */

  int GetUncorrectedEfficiency();
  /** gets the bbc efficinecy with the zvertex correction */
 private:
  int bbc_narrow_trig_bit_;
  int bbc_wide_trig_bit_;
  int zdc_wide_trig_bit_;
  // Triggerbits, hardcoded in constructor, can be overridden by calling member.
  // i.e Set___Trigger functions.

  std::map <int, int> bit_map_;
  /** breaks if you give negative input 
   * (which is incorrecte for integer valued return condition ).
   * Otherwise, will recursively multiply base times itself until
   * it is done exp times. */
  
  int ipow(int base, int exp);   
  int Derivative(TH1F*& h, TH1F*& d);
  /* Treats histogram, h, as a function. Returns derivative of h, stored
   * in a histogram, d. Remember, that root does not deal well with fitting
   * for histograms which have negative valued bins. */

  int FitVertexEfficiency(); 
  /** Generate first derivative from rebinned trigger acceptance plot. Since
   * we expect the trigger acceptance plot to be of the form of an error funtion
   * (i.e. the integral of a gaussian) we can fit the first derivative with two
   * gaussians, the peaks of which will represent the inflection points, which 
   * is the purpose of the fitting. This vertex range is then the range we will
   * use for the efficiency calculation. */



  std::string run_number_;
  std::string root_file_name_;
  std::string tree_name_;
  /* Vertex cut ranges */
  double z_vtx_cut_min_; /** the lower bound of the z-vertex cut */
  double z_vtx_cut_max_; /** the upper bound of the z-vertex cut */
  double z_vtx_full_max_;/** the upper bound of full z-vertex range */
  double z_vtx_full_min_;/** the lower bound of full z-vertex range */

  std::string BBCWideHistoName;
  std::string BBCNarrowHistoName;
  std::string BBCWideAndNarrowHistoName;
  std::string ZDCWideHistoName;
  std::string ZDCWideAndBBCWideHistoName;
  /** Histograms From Vernier Production */
  TFile* root_file_                ;
  TH1F* h_bbc_wide_z_            ;
  TH1F* h_bbc_narrow_z_          ;
  TH1F* h_bbc_wide_and_narrow_z_ ;
  TH1F* h_zdc_wide_z_            ;
  TH1F* h_zdc_and_bbc_z_         ;
  TH1F* trig_scaled_bits_        ;
  TH1F* trig_raw_bits_           ;
  TH1F* trig_live_bits_          ;
  /** fits are to abs_acceptance_first_derivative, generated 
   * from smoothed trigger acceptance trigger_acceptance_smoothed_ */
  TF1* left_gaus_;
  TF1* right_gaus_;

  /** Results From Efficiency Calculation */
  TH1F* vertex_efficiency_ ;
  TH1F* trigger_acceptance_;
  TH1F* trigger_acceptance_corrected_;
  TCanvas* c_trigger_acceptance_;
  TCanvas* c_vertex_efficiency_;
  TCanvas* c_corrected_Efficiency;
  TF1* vertex_fit_pol2_;
  TF1* vertex_fit_gaus_;
  TH1F* acceptance_first_derivative_;
  TH1F* abs_acceptance_first_derivative_;

  /** Display Results */
  TCanvas* c_BbcWide_zvtx              ;
  TCanvas* c_BbcNarrow_zvtx            ;
  TCanvas* c_ZdcWide_zvtx              ;
  TCanvas* c_BbcWideNarrowCoincidenceZ ;
  TCanvas* c_BbcZdcCoincidence_BbcZ    ;
  TCanvas* triggerAccCanvas            ;
  TCanvas* triggerVerCanvas            ;

  /** Store all histograms, canvases separately for drawing purposes */
  std::vector<TObject*> canvas_registry_;
  std::vector<TObject*> plot_registry_; // keep all objects to be drawn/saved in here.

};
#endif
