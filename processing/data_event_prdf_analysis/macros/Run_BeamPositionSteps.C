#include "FileManagement.h"
void Run_BeamPositionSteps(
    int run_index = 0,
    const std::string& vernierAnalysisLibrary = "libVernierTimeAnalysis.so") 
{
  gSystem->Load(vernierAnalysisLibrary.c_str());
  BeamPositionSteps bpm;
  bpm.Init(
    run_number[run_index],
    bpm_file[run_index],
    relative_time_step_boundaries[run_index]
     );
  bpm.Run();
  bpm.SaveEpochSteps("/direct/phenix+spin2/beaumim/vernierScans/prdf_analysis/prdf_tools/processing/time_compressed_processing/data");
  bpm.GetCentralStepTime(3);
  bpm.MakeFigures("/direct/phenix+spin2/beaumim/vernierScans/prdf_analysis/prdf_tools/processing/time_compressed_processing/plots");
}
