int run_packet_dump(
    const std::string& PRDF_name = "/gpfs02/phenix/spin3/beaumim/Run12VernierPRDF/files/VSCANDATA_P00-0000364636-0005.PRDFF",
    const std::string& out_file_name = "testout.txt",
    const std::string& lib = "/direct/phenix+spin2/beaumim/vernierScans/prdf_analysis/prdf_tools/scaler_event_prdf_scalers/libclock.so"
)
{
  // Options
  // pstart() - call this to run over the whole PRDFF, you can still enter commands from root
  // prun(0) - call this to run without option to stop
  // pstop() - stop interating over PRDFF. See packetmonitor.h for more details ($ONLINE_MAIN/includes/packetmonitor.h)
  gSystem->Load(lib.c_str());


  set_output_file(out_file_name);
  pfileopen(PRDF_name.c_str());
  prun(0);
  finish();
  std::cout << "All done." << std::endl;
  return 0;
}
