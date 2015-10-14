#include "../../../FileManagement.h"

int Run_HourglassData(
  int run_index = 0,
  const std::string lib = "libVernierHourglass.so"
) {
  gSystem->Load(lib.c_str());
  HourglassData hd;
  hd.Init(
    run_number[run_index],
    merged_prdf_dst[run_index],
    epoch_step_boundaries[run_index],
    bpm_file[run_index],
    relative_time_step_boundaries[run_index]
  );
  hd.Run();
  hd.SaveFigures("/direct/phenix+spin2/beaumim/vernierScans/plots");
  return 0;
}
