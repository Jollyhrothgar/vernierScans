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
      const std::string& epoch_steps_file,
      const std::string& bpm_data_file_name, 
      const std::string& bpm_steps_file_name
      );
  int Run();
  int SaveFigures(const std::string& output_dir); 
 private:
  std::string run_number_;
  int LoadSteps();
  int InitHistograms();
  int CreatePlots(); /** Create some visually pleasing plots */

  BeamPositionSteps bpm_;
  std::vector<TObject*> plot_registry_;
  std::vector<std::pair<long int, long int> > step_boundaries_;
  std::string epoch_steps_file_name_;
  std::string scalers_file_name_;
  std::map<int,TH1F*> bbc_z_vtx_; /** z-vertex distribution for each step */
  std::map<int,TH1F*> zdc_z_vtx_; /** z-vertex distribution for each step */
  std::vector<double> h_steps_;
  std::vector<double> v_steps_;
  std::vector< std::pair< double, bool > > displacement_; /** maps total beam displacement to bool 1 = mostly horizontal step, 0 = mostly vertical step */
};

#endif
