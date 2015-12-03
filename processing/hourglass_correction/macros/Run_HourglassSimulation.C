#include "../../../FileManagement.h"

int Run_HourglassSimulation(
  const std::string& file_stub,
  const int simulation_mode = 0,
	const int zprofile_mode = 0,
  const std::string lib = "libVernierHourglass.so"
) {
  int run_index = 0;
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  std::string config_file = file_stub + ".conf";
  std::string compare_file = file_stub + "_compare.root";
  sim.InitFromConfig(config_file);
  sim.OverrideSaveFile(file_stub); // Use for batch jobs where complex file name would be confusing.
  sim.LoadZProfile(
    z_profile_density_blue[run_index],
    z_profile_density_yellow[run_index]
  );

  switch (simulation_mode) {
    case 0:
      sim.Run(zprofile_mode);
      sim.Compare(hourglass_data_file[run_index]);
      sim.SaveFigures(plots_dir);
      break;
    case 1:
      sim.RunRootFinder(zprofile_mode,hourglass_data_file[run_index]);
      break;
    default:
      std::cout << "no adequate running mode supplied." << std::endl;
      std::cout << " choose 1 for RootFinder mode" << std::endl;
      std::cout << " choose 0 for Normal Running mode " << std::endl;
  } 
  return 0;
}
