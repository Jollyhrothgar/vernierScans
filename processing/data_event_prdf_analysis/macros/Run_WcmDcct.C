#include "../../../FileManagement.h"

void Run_WcmDcct(
    int run_index = 0,
    const std::string& lib = "libVernierTimeAnalysis.so"
    )
{
  gSystem->Load(lib.c_str());
  WcmDcctManager wdm;
  wdm.Init(
      run_number[run_index],
      wcm_blue_file[run_index],
      wcm_yellow_file[run_index],
      dcct_file[run_index]
      );
  BeamPositionSteps bpm;
  bpm.Init(
    run_number[run_index],
    bpm_file[run_index],
    relative_time_step_boundaries[run_index]
     );
  bpm.Run();
  float scan_seconds = (float)bpm.GetScanEndTime() - (float)bpm.GetScanStartTime();
  cout << "The whole scanning portion of the run took: " 
      << scan_seconds << " seconds (" << scan_seconds / 3600.00 << ") hours" << std::endl;
  wdm.SetTimeOffset(bpm.GetTimeOffset());
  wdm.SetScanRange(bpm.GetScanStartTime(), bpm.GetScanEndTime());
  wdm.Run();
  double lumi_pct_loss = wdm.GetLuminosityLoss();
  std::cout << "Lost " << lumi_pct_loss*100.0 << " percent luminosity!" << std::endl;
  wdm.SaveFigures(plots_dir);
  wdm.SaveBeamPopulations(summary_dir);
  //wdm.PrintDataTimeIndex(1330279088);
  //wdm.PrintDataTimeIndex(1330279088);
  //wdm.PrintBlueWcmTotal();
  //wdm.PrintYellowWcmTotal();
  //wdm.PrintYellowWcmTotal();
}                                                                                                    
