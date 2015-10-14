#include "../../../FileManagement.h"
int Run_PRDFBBCRate(
  int run_index = 0,
  const std::string& lib = "libPRDFAnalysis.so"
) {
  gSystem->Load(lib.c_str());
  PRDFBBCRate prdf_rate;
  prdf_rate.Init(
    run_number[run_index],
    scaler_ppg_file[run_index]
    );
  prdf_rate.Run();
  prdf_rate.SaveFigures("/direct/phenix+spin2/beaumim/vernierScans/prdf_analysis/prdf_tools/processing/scaler_event_processing/plots");
  return 0;
}
