#include "FileManagement.h"

int Run_HourglassSimulation(
  int run_index = 0,
  const std::string lib = "libVernierHourglass.so"
) {
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  sim.Init();
  sim.Run();
  sim.SaveFigures("/direct/phenix+spin2/beaumim/vernierScans/plots");
  return 0;
}
	
