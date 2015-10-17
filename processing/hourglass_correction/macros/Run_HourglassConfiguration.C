#include "../../../FileManagement.h"

int Run_HourglassConfiguration(
  int run_index = 0,
  const std::string lib = "libVernierHourglass.so"
) {
  gSystem->Load(lib.c_str());
  HourglassConfiguration con;
  // Create All Config Files
  con.ModifyConfigParameter("RUN_NUMBER","359711");
  con.SetDefaultValues();
  con.ShowConfigFile();
  con.SaveConfigFile("/direct/phenix+spin2/beaumim/vernierScans/data/run_12/simulation_config");
  return 0;
}
