#include "../../../FileManagement.h"

int Run_HourglassSimulation(
  int run_index = 0,
  const std::string lib = "libVernierHourglass.so"
) {
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  //sim.InitDefault();
  sim.InitFromConfig("/direct/phenix+spin2/beaumim/vernierScans/processing/hourglass_correction/macros/test_config.conf");
  sim.OverrideSaveFile("./test_simulation");
  sim.Run();
  sim.SaveFigures(plots_dir);
  sim.Compare(hourglass_data_file[run_index],"./test_compare.root");
  return 0;
}
