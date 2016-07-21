#ifndef __BEAM_POSITION_STEPS__
#define __BEAM_POSITION_STEPS__

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>
#include <cmath>

#include "TH1F.h"
#include "TGraph.h"

#include "BeamSeparationData.h"

/* Class currently configured to treat the first time entry
 * in the data file as "0" seconds. This was done so as to be
 * able to identify relative times by eye on data plotted
 * against time (since epoch time is such a large number, it cannot
 * be displayed on the x-axis conveniently). 
 * Beam steps are therefore encoded in this relative time, and 
 * must be modified with the correct offset to be used with data 
 * using the fill data. The offset is printed when this class is used. */
class BeamPositionSteps
{
  public:
   /*** DST ANALYSIS ***/
   BeamPositionSteps();
   ~BeamPositionSteps();
   int Init(const std::string& runNumber, const std::string& bpmDataFileName, const std::string& stepsFileName);
   int Run();
  private:
   int DrawScan();
  public:
   float GetHStep(int i);
   float GetVStep(int i);
   float GetHStepErr(int i);
   float GetVStepErr(int i);
   unsigned int GetNumberOfSteps(); 
   /** returns size of the vernierScanBPMSteps object */

   time_t GetTimeOffset(); 
   /** returns the time corresponding to the first entry of the data file. Use
    * to sync with WcmDcct Data */

   time_t GetScanStartTime();
   time_t GetScanEndTime();
   time_t GetCentralStepTime(int step_i ); 
   /** Returns the epoch time closest to the center of the step# */

   // Estimate BPM accuracy by looking at RMS of fill, corrected for beam
   // movement
   int SaveBeamPositionData(const std::string& out_file_name);

   /*** PRDF ANALYSIS ***/ 
   int LookupStep(double lookup_time); 
   /** takes an epoch time stamp and determines which scan step it belongs to.
    * Double is used for half-second precision in cases where we need to
    * average betwee time stamps. */

  int SaveEpochSteps(const std::string& out_file_directory);
  /** takes the beam steps file that was loaded in LoadSteps and dumps epoch
   * time to a different file. This is useful so that we don't have to drag
   * around some offset time everywhere we want to compare step boundaries.
   * Step boundaries based on BPM data will be most conservative, since BPM
   * data rate is once every four seconds, as opposed to BBC data dumped from
   * PRDFs, which is every second, or faster if event-sequence is used as a
   * time proxy.
   */
   
   double bpm_global_rms_;
   double bpm_global_average_rms_;
   std::map<long int, BeamSeparationData> beam_separation_data_;
   std::map<long int, BeamSeparationData> GetBeamSeparationData();
   /** beam_separation_data_ and GetBeamSeparationData : organize bpm data
    * relevant to beam width into one structure All uncertainties are
    * determined by the average RMS of the beam for each step, but the actual
    * data in the beam_separation_data_ object is organized in time.
    * BeamSeparationData fills bpm_rms into BeamSeparationObject stored in
    * beam_separation_data_ and returns the whole beam_separation_data_ object.
    * Data is binned in time */

   int BeamPositionLookup(double time_lookup, double rejection_threshold, double& x_sep, double& y_sep); 
   /** Intripolate the nearest beam position separation based on a time-stamp
    * lookup Return 0 if successful, 1 if not. This is used for raw bbc data
    * dumped to a text file. x_sep and y_sep are modified, so you have to
    * declare these variables ahead of time. rejection_threshold sets a
    * threshold difference for the lookup - if the found time is
    * rejection_threshold seconds away from the lookup, then the
    * BeamPositionLookup will retun "1". When points match exactly, they are
    * returned, otherwise we do a linear intripolation. */
   
   bool IsFirstScan(time_t time_index);
   /** These functions will determine, from the steps files, whether or not a
    * time index is in the first scan sweep */

   bool IsSecondScan(time_t time_index);
   /** These functions will determine, from the steps files, whether or not a
    * time index is in the second scan sweep */

  private:
   int LoadBpmData();
   int MakeSteps(); /** Generate steps from data */
   bool error_state;
   
   // private nested struct to hide implementation of data organization to keep
   // namespace clear and to focus user on public interface of this class.
   //
   // bl[xy]Avg and ye[xy]Avg are the average of the BPMs' measurement of x and
   // y, with hSeparationAvg and vSeparationAvg represeting the beam separation
   // at the PHENIX IR calculated from the averages
   //
   // {bl,ye}{x,y}IR are the x and y positions at the IR calculated via the
   // linear extrapolation of the beams between the stations, with hSeparation
   // and vSeparation representing the separations calculated from these values
   struct BeamPosition
   {
     float blx7;
     float bly7;
     float yex7;
     float yey7;
     float blx8;
     float bly8;
     float yex8;
     float yey8;
     float blxIR;
     float blyIR;
     float yexIR;
     float yeyIR;
     float blxAvg;
     float blyAvg;
     float yexAvg;
     float yeyAvg;
     time_t time;
     float hSeparation;
     float vSeparation;
     float hSeparationAvg;
     float vSeparationAvg;

     int Init() {
       bly7  = -999;
       yex7  = -999;
       yey7  = -999;
       blx8  = -999;
       bly8  = -999;
       yex8  = -999;
       yey8  = -999;
       blxIR = -999;
       blyIR = -999;
       yexIR = -999;
       yeyIR = -999;
       blxAvg = -999;
       blyAvg = -999;
       yexAvg = -999;
       yeyAvg = -999;
       hSeparationAvg = -999;
       vSeparationAvg = -999;
       time  = -999;
       hSeparation = -999;
       vSeparation = -999;
       return 0;
     }

     void CalculateIRPositions() {
       // hSeparation is mathematically equivalent to hSeparationAvg!!!
       // likewise for the rest.
       //
       // hSeparation == blx7+0.5*(blx8 - blx7) - (yex7+0.5*(yex8 - yex7))
       //      == blx7 + 0.5*blx8 - 0.5*blx7 - yex7 - 0.5*yex8 + 0.5*yex7
       //      == 0.5*blx7 + 0.5*blx8 - (0.5*yex7+0.5*yex8)
       //      But 0.5*blx7+0.5*blx8 = 0.5*(blx7+blx8) == blxAvg !!! And so on.
       blxIR = blx7+0.5*(blx8 - blx7);
       blyIR = bly7+0.5*(bly8 - bly7);
       yexIR = yex7+0.5*(yex8 - yex7);
       yeyIR = yey7+0.5*(yey8 - yey7);
       blxAvg = 0.5*(blx8+blx7);
       blyAvg = 0.5*(bly8+bly7);
       yexAvg = 0.5*(yex8+yex7);
       yeyAvg = 0.5*(yey8+yey7);
       hSeparation = blxIR - yexIR;
       vSeparation = blyIR - yeyIR;
       hSeparationAvg = blxAvg - yexAvg;
       vSeparationAvg = blyAvg - yeyAvg;
     }
     void Print() {
       std::cout << " " << time 
           << " " << blx7 << " " << bly7 << " " << yex7 << " " << yey7 
           << " " << blx8 << " " << bly8 << " " << yex8 << " " << yey8 
           << " " << blxIR << " " << blyIR << " " << yexIR << " " << yeyIR 
           << " " << hSeparation << " " << vSeparation 
           << " " << blxAvg << " " << blyAvg << " " << yexAvg << " " << yeyAvg 
           << " " << hSeparationAvg << " " << vSeparationAvg << std::endl; 
     }
   };

   // private nested class to hide organization/implementation of BPM step
   // data.
   // This struct holds the real vernier scan BPM separation
   struct BeamStep
   {
     BeamStep() {
       bpmUncertainty = 10.;
     }
     // directly assign these variables.
     
     // Calculated from linear extrapolation
     float xSeparation;
     float ySeparation;
     float xSeparationRMS;
     float ySeparationRMS;

     // Calculated from average of station 7 and 8
     float xSeparationAvg;
     float ySeparationAvg;
     float xSeparationAvgRMS;
     float ySeparationAvgRMS;

     float bpmUncertainty; // microns, PHENIX technical note: tn406.0
     float GetXSeparation() {
       return xSeparation;
     }
     float GetYSeparation() {
       return ySeparation;
     }
     float GetXSeparationError() {
       return pow(pow(xSeparationRMS,2.0) + pow(bpmUncertainty,2.0),0.5);
     }
     float GetYSeparationError() {
       return pow(pow(ySeparationRMS,2.0) + pow(bpmUncertainty,2.0),0.5);
     }
   };

   std::vector<BeamStep> vernierScanBPMSteps;
   std::vector<BeamStep> beamSteps;
   int LoadStepsFile();

   std::string bpmDataFileName;
   std::string stepsFileName;
   std::string runNumber;
   // maps c time stamp to data set
   std::map<time_t, BeamPosition> bpmData;
   std::vector<std::pair< int, int > > steps; 
   /* holds the beginning and end of each scan step. Add GetTimeOffset() to get
    * the epoch time of the step. steps stores relative time to start of BPM
    * data file for the run (not the whole fill). This complicated crap is
    * because otherwise, its hard to tell where in the scan we are by eye since
    * we'd have an epoch stamp, as opposed to a simple number between 0 and the
    * number of seconds in the bpm data file for the run 
    *
    * Also, we assume an even number of steps */
  public:
   TGraph* beam_separation_v;
   TGraph* beam_separation_h; 
  private:
   std::vector<TObject*> plot_registry_;
  public:
   int MakeFigures(const std::string& figure_output_dir);
};

#endif
