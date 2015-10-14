#ifndef __MERGED_PRDF_DST__
#define __MERGED_PRDF_DST__
struct MergedPRDFDST{
  int triglive;
  int trigscaled;
  int trigraw;
  time_t timestamp;
  float bbc_z;
  float zdc_z;
  int bunch;
  int gl1p_bbc;
  int atp_number;
  int gl1p_clock;
  int gl1p_zdc_wide;
  int gl1p_zdc_narrow;

  void Reset(){
    triglive = 0;
    trigscaled = 0;
    trigraw = 0;
    timestamp = 0;
    bbc_z = 0;
    zdc_z = 0;
    bunch = 0;
    gl1p_bbc = 0;
    gl1p_clock = 0;
    gl1p_zdc_wide = 0;
    gl1p_zdc_narrow = 0;
  }
};
#endif
