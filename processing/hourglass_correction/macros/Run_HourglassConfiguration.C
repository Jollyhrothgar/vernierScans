#include "../../../FileManagement.h"

int Run_HourglassConfiguration(
) {
  gSystem->Load("libVernierHourglass.so");
  gSystem->Load("libVernierTimeAnalysis.so");
  HourglassConfiguration con;
  // Create All Config Files
  for(int run_index = 0; run_index < NUMBER_OF_RUNS; run_index++) {
    hd.Run();
  }
   
  con.ModifyConfigParameter("RUN_NUMBER","359711");
  con.SetDefaultValues();
  con.ShowConfigFile();
  con.SaveConfigFile("/direct/phenix+spin2/beaumim/vernierScans/data/run_12/simulation_config");
  return 0;
}
