#include "../../../FileManagement.h"
void Run_BeamPositionSteps(
  int run_index = 0,
  const std::string& vernierAnalysisLibrary = "libPRDFAnalysis.so") 
{
  gSystem->Load(vernierAnalysisLibrary.c_str());
  BeamPositionSteps bpm;
  bpm.Init(
    run_number[run_index],
    bpm_file[run_index],
    steps_boundaries[run_index]
  );
  bpm.Run();
  bpm.GetCentralStepTime(3);
}
