#include "../../../FileManagement.h"

int Run_HourglassSimulation(
  const std::string& file_stub,
  const int simulation_mode = 0,
  const int zprofile_mode = 0,
  const std::string& plot_output = plots_dir,
  const std::string lib = "libVernierHourglass.so"
) {
  int run_index = 0;
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  std::string config_file = file_stub + ".conf";
  sim.SetSaveFileStub(file_stub);
  sim.Init(
      config_file,
      hourglass_data_file[run_index],
      z_profile_density_blue[run_index],
      z_profile_density_yellow[run_index]
      );
  sim.SetSaveDirectory(plot_output);
  sim.SetSaveAll();
  switch (simulation_mode) {
    case 0:
      sim.Run(zprofile_mode);
      sim.Compare();
      break;
    case 1:
      sim.RunRootFinder(zprofile_mode,hourglass_data_file[run_index]);
      break;
    default:
      std::cout << "no adequate running mode supplied." << std::endl;
      std::cout << " choose 1 for RootFinder mode" << std::endl;
      std::cout << " choose 0 for Normal Running mode " << std::endl;
  } 
  sim.SaveFigures();
  return 0;
}
