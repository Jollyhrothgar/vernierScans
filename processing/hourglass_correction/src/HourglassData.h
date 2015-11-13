#ifndef __HOURGLASS_DATA_H__
#define __HOURGLASS_DATA_H__

// Header File For Linked Libraries
#include "/direct/phenix+spin2/beaumim/vernierScans/processing/dst_analysis/src/BeamPositionSteps.h"

// STL
#include <string>
#include <vector>
#include <map>
#include <algorithm>

// ROOT Data
#include "TH1F.h"
#include "TObject.h"

class HourglassData {
 public:
  HourglassData();
  ~HourglassData();
  int Init(
      const std::string& run_number,
      const std::string& scalers_file,
      const std::string& epoch_step_boundaries,
      const std::string& bpm_data_file_name, 
      const std::string& relative_step_boundaries,
      const std::string& bpm_planned_steps_file_name
      );
  int Run();
  int SaveFigures(const std::string& output_dir); 
  int SaveSimulationConfigData(const std::string& output_h);
 private:
  std::string run_number_;

  // This loads the epoch times which bound each scan step. We use this after
  // merging the bbc dst data with the prdf time stamped data.
  int LoadEpochStepBoundaries();
  int LoadPlannedSteps();
  int InitHistograms(); 

  // based on average value of histograms, show the average offset of the BBC
  // and ZDC z-vertex 
  int ShowOffsets(); 
 
  BeamPositionSteps bpm_;
  std::vector<TObject*> plot_registry_;
  std::vector<std::pair<long int, long int> > step_boundaries_;
  
  // CAD planned total beam displacements.
  std::vector<float> planned_steps_; 
  std::string epoch_step_boundaries_file_name_;
  std::string bpm_planned_steps_file_name_;
  std::string scalers_file_name_;

  // z-vertex distribution for each step 
  std::map<int,TH1F*> bbc_z_vtx_; 

  // z-vertex distribution for each step 
  std::map<int,TH1F*> zdc_z_vtx_; 
  std::vector<double> h_steps_;
  std::vector<double> v_steps_;

  // The offset between the bbc zvertex reckoning and zdc z-vertex reckoning. If
  // many offsets are availalble, the average is used.
  float bbc_zdc_offset_;

  //  maps total beam displacement to bool 1 = mostly horizontal step, 0 =
  //  mostly vertical step
  std::vector< std::pair< double, bool > > displacement_;
};

#endif
