#include "../../../FileManagement.h"

int Run_BeamWidthTime(
  int run_index = 0,
  int bunch_index = -1,
  const std::string& vernierAnalysisLibrary = "libVernierTimeAnalysis.so") 
{
  gSystem->Load(vernierAnalysisLibrary.c_str());
  char run [256];
  char file[256];
  if( bunch_index >= 0 && bunch_index < 120) {
    sprintf(run,"%s_bunch_%d",run_number[run_index].c_str(),bunch_index);
    printf("%s\n",run);
    sprintf(file,"%s_%d.txt",scaler_file_bunch_stub[run_index].c_str(),bunch_index);
    printf("%s\n",file);
  } else {
    sprintf(run,"%s",run_number[run_index].c_str());
    printf("%s\n",run);
    sprintf(file,"%s",scaler_file[run_index].c_str());
    printf("%s\n",file);
  }
  BeamWidthTime bws;
  bws.Init(
      run,
      bpm_file[run_index],
      relative_time_step_boundaries[run_index],
      epoch_step_boundaries[run_index], 
      file,
      planned_beam_steps[run_index]
      );
  bws.Run();
  bws.SaveFigures("/direct/phenix+spin2/beaumim/vernierScans/plots");
  bws.SaveBeamWidthData("/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data");
  bws.MakePlannedStepsTable("/direct/phenix+WWW/p/draft/beaumim/Vernier_Analysis/run_12_vernier_analysis_note/vernier/BeamPositionMonitoring/tables","plannedsteps");
  return 0;
}
