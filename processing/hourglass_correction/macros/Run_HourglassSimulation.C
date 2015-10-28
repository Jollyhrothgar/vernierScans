#include "../../../FileManagement.h"

int Run_HourglassSimulation(
  const std::string& config_file,
  const std::string& save_file_stub,
  const std::string lib = "libVernierHourglass.so"
) {
  int run_index = 0;
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  //sim.InitDefault();

  std::string compare_file = save_file_stub + "_compare.root";
  //sim.InitFromConfig(simulation_config_dir+"/359711_step_00.conf");
  sim.InitFromConfig(config_file);
  sim.OverrideSaveFile(save_file_stub); // Use for batch jobs where complex file name would be confusing.
  sim.Run();
  sim.SaveFigures(plots_dir);
  sim.Compare(hourglass_data_file[run_index],compare_file);
  return 0;
}
