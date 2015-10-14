#include "FileManagement.h"
int Run_BeamWidthScalers(
  int run_index = 0,
  const std::string& vernierAnalysisLibrary = "libPRDFAnalysis.so") 
{
  gSystem->Load(vernierAnalysisLibrary.c_str());
  BeamWidthScalers bws;
  bws.Init(
    run_number[run_index],
    bpm_file[run_index],
    steps_boundaries[run_index],
    scaler_ppg_file[run_index] 
    );
  bws.Run();
  bws.Draw();
  return 0;
}
