#include "../../../FileManagement.h"

int Run_HourglassSimulation(
  const std::string& file_stub,
  const std::string lib = "libVernierHourglass.so"
) {
  int run_index = 0;
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  //sim.InitDefault();
  std::string config_file = file_stub + ".conf";
  std::string compare_file = file_stub + "_compare.root";
  sim.InitFromConfig(config_file);
  sim.OverrideSaveFile(file_stub); // Use for batch jobs where complex file name would be confusing.
  sim.Run();
  sim.Compare(hourglass_data_file[run_index]);
  sim.SaveFigures(plots_dir);
  return 0;
}
