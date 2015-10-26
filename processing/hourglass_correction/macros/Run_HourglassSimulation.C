#include "../../../FileManagement.h"

int Run_HourglassSimulation(
  int run_index = 0,
  const std::string lib = "libVernierHourglass.so"
) {
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  //sim.InitDefault();

  sim.InitFromConfig(simulation_config_dir+"/359711_step_00.conf");
  sim.OverrideSaveFile("./359711_step_00_vertex");
  sim.Run();
  sim.SaveFigures(plots_dir);
  sim.Compare(hourglass_data_file[run_index],"./359711_step_00.root");
  return 0;
}
