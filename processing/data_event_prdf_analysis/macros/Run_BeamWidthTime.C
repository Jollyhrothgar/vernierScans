#include "../../../FileManagement.h"

int Run_BeamWidthTime(
  int run_index = 0,
  int bunch_index = -1,
  const std::string& vernierAnalysisLibrary = "libVernierTimeAnalysis.so") 
{
  gSystem->Load(vernierAnalysisLibrary.c_str());
  std::stringstream run;
  std::stringstream file;
  if( bunch_index >= 0 && bunch_index < 120) {
    run.str("");
    file.str("");
    run << run_number[run_index] << "_bunch_" << bunch_index;
    file << scaler_file_bunch_stub[run_index] << bunch_index << ".txt";
  } else {
    run.str("");
    file.str("");
    run << run_number[run_index];
    file << scaler_file[run_index]; 
  }
  BeamWidthTime bws;
  bws.Init(
      run.str(),
      bpm_file[run_index],
      relative_time_step_boundaries[run_index],
      epoch_step_boundaries[run_index], 
      file.str(),
      planned_beam_steps[run_index]
      );
  bws.Run();
  bws.SaveFigures("/direct/phenix+spin2/beaumim/vernierScans/plots");
  bws.SaveBeamWidthData("/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data");
  bws.MakePlannedStepsTable("/direct/phenix+WWW/p/draft/beaumim/Vernier_Analysis/run_12_vernier_analysis_note/vernier/BeamPositionMonitoring/tables","plannedsteps");
  return 0;
}
