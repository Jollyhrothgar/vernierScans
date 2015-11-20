#ifndef __FILE_MANAGEMENT__
#define __FILE_MANAGEMENT__

// Note, dealing with CINT's limitations means that we cannot map these values
// to some kind of nTuple or a std::map, so its imparative that we maintain the proper 
// order for every array here. Note that, in all cases, the run order is maintained.
// There is nothing special about the run order, except that the same index in all arrays
// access data for the same run.

const int NUMBER_OF_RUNS = 7;


// Global File Output
const std::string analysis_root = "/direct/phenix+spin2/beaumim/vernierScans";
const std::string simulation_config_dir = "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/simulation_config";
const std::string plots_dir = "/direct/phenix+spin2/beaumim/vernierScans/plots";
const std::string summary_dir = "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data";

std::string run_number[11] = {
  "359711",
  "360879",
  "362492",
  "364636",
  "365866",
  "366605",
  "367138",
  "431624",
  "431723",
  "431857",
  "431962"
};

std::string bpm_fill_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bpm_fill_16444.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bpm_fill_16470.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bpm_fill_16514.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bpm_fill_16587.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bpm_fill_16625.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bpm_fill_16655.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bpm_fill_16671.dat"
};

std::string bpm_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/359711_BPM.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/360879_BPM.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/362492_BPM.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/364636_BPM.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/365866_BPM.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/366605_BPM.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/367138_BPM.dat"
};
std::string wcm_yellow_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/359711_WCM_16444_Yellow.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/360879_WCM_16470_Yellow.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/362492_WCM_16514_Yellow.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/364636_WCM_16587_Yellow.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/365866_WCM_16625_Yellow.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/366605_WCM_16655_Yellow.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/367138_WCM_16671_Yellow.dat"
};
std::string wcm_blue_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/359711_WCM_16444_Blue.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/360879_WCM_16470_Blue.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/362492_WCM_16514_Blue.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/364636_WCM_16587_Blue.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/365866_WCM_16625_Blue.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/366605_WCM_16655_Blue.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/367138_WCM_16671_Blue.dat"
};
std::string relative_time_step_boundaries[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/359711_steps.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/360879_steps.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/362492_steps.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/364636_steps.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/365866_steps.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/366605_steps.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/367138_steps.dat"
};
std::string bbc_bin_step_boundaries[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/359711_bbcrate_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/360879_bbcrate_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/362492_bbcrate_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/364636_bbcrate_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/365866_bbcrate_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/366605_bbcrate_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/367138_bbcrate_steps.txt"
};
std::string epoch_step_boundaries[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/359711_epoch_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/360879_epoch_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/362492_epoch_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/364636_epoch_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/365866_epoch_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/366605_epoch_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/367138_epoch_steps.txt"
};
std::string dcct_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/359711_DCCT_16444.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/360879_DCCT_16470.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/362492_DCCT_16514.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/364636_DCCT_16587.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/365866_DCCT_16625.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/366605_DCCT_16655.dat",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/367138_DCCT_16671.dat"
};
std::string scaler_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/359711_bunch_integrated.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/360879_bunch_integrated.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/362492_bunch_integrated.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/364636_bunch_integrated.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/365866_bunch_integrated.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/366605_bunch_integrated.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/367138_bunch_integrated.txt"
};
std::string scaler_file_bunch_stub[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/359711_bunch_",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/360879_bunch_",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/362492_bunch_",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/364636_bunch_",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/365866_bunch_",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/366605_bunch_",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_data/367138_bunch_"
};
std::string scaler_ppg_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/scaler_data/ppg_scalers_359711.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/scaler_data/ppg_scalers_360879.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/scaler_data/ppg_scalers_362492.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/scaler_data/ppg_scalers_364636.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/scaler_data/ppg_scalers_365866.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/scaler_data/ppg_scalers_366605.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/scaler_data/ppg_scalers_367138.txt"
};
std::string merged_prdf_dst[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_dst_merged/359711_dst_prdf.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_dst_merged/360879_dst_prdf.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_dst_merged/362492_dst_prdf.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_dst_merged/364636_dst_prdf.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_dst_merged/365866_dst_prdf.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_dst_merged/366605_dst_prdf.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/prdf_dst_merged/367138_dst_prdf.txt"
};
std::string reduced_dst_file[11] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/359711_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/360879_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/362492_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/364636_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/365866_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/366605_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/367138_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_15/431624/431624_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_15/431723/431723_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_15/431857/431857_reduced.root",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_15/431962/431962_reduced.root"
};
std::string planned_beam_steps[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/359711_planned_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/360879_planned_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/362492_planned_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/364636_planned_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/365866_planned_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/366605_planned_steps.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/367138_planned_steps.txt"
};

std::string beam_width_summary_data[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_BeamWidthData.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_BeamWidthData.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_BeamWidthData.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_BeamWidthData.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_BeamWidthData.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_BeamWidthData.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_BeamWidthData.txt"
};

std::string zdc_counts_per_step[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_ZDCCountsPerStep.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_ZDCCountsPerStep.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_ZDCCountsPerStep.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_ZDCCountsPerStep.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_ZDCCountsPerStep.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_ZDCCountsPerStep.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_ZDCCountsPerStep.txt"
};

std::string x_offset_sim_config[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_XOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_XOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_XOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_XOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_XOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_XOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_XOffsets.txt"
};

std::string y_offset_sim_config[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_YOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_YOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_YOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_YOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_YOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_YOffsets.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_YOffsets.txt"
};

std::string v_width_sim_config[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_vWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_vWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_vWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_vWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_vWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_vWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_vWidth.txt"
};

std::string h_width_sim_config[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_hWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_hWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_hWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_hWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_hWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_hWidth.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_hWidth.txt"
};

std::string beam_population_sim_config[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_WCMDCCT_BeamPopulation.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_WCMDCCT_BeamPopulation.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_WCMDCCT_BeamPopulation.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_WCMDCCT_BeamPopulation.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_WCMDCCT_BeamPopulation.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_WCMDCCT_BeamPopulation.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_WCMDCCT_BeamPopulation.txt"
};

std::string zdc_bbc_offset_sim_config[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_BBCZDCOffset.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_BBCZDCOffset.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_BBCZDCOffset.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_BBCZDCOffset.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_BBCZDCOffset.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_BBCZDCOffset.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_BBCZDCOffset.txt"
};

std::string zdc_zvtx_histo_name_sim_config[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/359711_ZDCVertexPlots.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/360879_ZDCVertexPlots.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/362492_ZDCVertexPlots.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/364636_ZDCVertexPlots.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/365866_ZDCVertexPlots.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/366605_ZDCVertexPlots.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/summary_data/367138_ZDCVertexPlots.txt"
};

std::string hourglass_data_file[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/plots/359711_HourglassData.root",
  "/direct/phenix+spin2/beaumim/vernierScans/plots/360879_HourglassData.root",
  "/direct/phenix+spin2/beaumim/vernierScans/plots/362492_HourglassData.root",
  "/direct/phenix+spin2/beaumim/vernierScans/plots/364636_HourglassData.root",
  "/direct/phenix+spin2/beaumim/vernierScans/plots/365866_HourglassData.root",
  "/direct/phenix+spin2/beaumim/vernierScans/plots/366605_HourglassData.root",
  "/direct/phenix+spin2/beaumim/vernierScans/plots/367138_HourglassData.root"
};

std::string z_profile_density_blue[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bwcm_zprofile_359711_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bwcm_zprofile_360879_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bwcm_zprofile_362492_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bwcm_zprofile_364636_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bwcm_zprofile_365866_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bwcm_zprofile_366605_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/bwcm_zprofile_367138_density.txt"
};

std::string z_profile_density_yellow[7] = {
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/ywcm_zprofile_359711_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/ywcm_zprofile_360879_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/ywcm_zprofile_362492_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/ywcm_zprofile_364636_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/ywcm_zprofile_365866_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/ywcm_zprofile_366605_density.txt",
  "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/cad_data/ywcm_zprofile_367138_density.txt"
};
#endif
