#include <iostream>
#include <fstream>

#include <pmonitor.h>
#include <packet_gl1.h>

#include "EventTypes.h"
#include "run_scalers.h"

/* Packet ID, annotations
 * 14001 - crossing information, beam polarization, spin patterns
 * 14009 - this seems to be the clock counter for the current, and next three events.
 * 14008 - this seems to be the gl1-1p scaler packet, there are two rows with four entries each
 * 14011 - this seems to have the crossing number, and additional information. */
const int gl1_packet_number     = 14001;  
const int bbc_ll1_packet_number = 14002;
const int gl1p1_packet_number   = 14008;
const int gl1p2_packet_number   = 14009;
const int gl1psum_packet_number = 14011;


// This is the GL1P board configuration for Run 12 Vernier Analysis
const int gl1p_scaler_id_bbc        = 0;
const int gl1p_scaler_id_clock      = 1;
const int gl1p_scaler_id_zdc_wide   = 2;
const int gl1p_scaler_id_zdc_narrow = 3;

static std::ofstream out_file;
int init_done = 0;

int set_output_file(const std::string& out_file_name) {
  out_file.open(out_file_name.c_str());
  return 0;
}

int finish() {
  out_file.close();
  return 0;
}

int pinit()
{
  if (init_done) return 1;
  init_done = 1;
  return 0;
}

int process_event (Event * e)
{
  // Get the ATP Number which processed/compressed the event
  int atp_number = e->getFrameEntry("FRAMESOURCEID");
  long long unsigned int time = e->getTime();
  /* The time-stamp needs to be fixed using the ATP calibration, but
   * here, it will just get spat-out into the output stream. The file
   * will then by synchronized, sorted, cleaned, in no particular
   * order. */

  
  if(e->getEvtType() == DATAEVENT){
    Packet* p_gl1p     = e->getPacket( gl1p1_packet_number   );
    Packet* p_gl1      = e->getPacket( gl1_packet_number     );
    Packet* p_bbcll1   = e->getPacket( bbc_ll1_packet_number );
    
    // Event number / event sequence / etc
    long unsigned int event_number = e->getEvtSequence();

    // Crossing number 
    int gl1_crossing_id = p_gl1 -> iValue(0,"CROSSCTR");
    gl1_crossing_id = (gl1_crossing_id+5)%120; // correct for standard crossing shift
    // IDENTICAL: int gl1_crossing_id2 = p_gl1p -> iValue(0,"CLOCK");
    
    // these are gl1-1p scalres, so they do not need to be corrected for overflow.
    long long unsigned int gl1p_scaler_bbc        = p_gl1p->iValue(gl1p_scaler_id_bbc        ,  "SCALER");
    long long unsigned int gl1p_scaler_clock      = p_gl1p->iValue(gl1p_scaler_id_clock      ,  "SCALER");
    long long unsigned int gl1p_scaler_zdc_wide   = p_gl1p->iValue(gl1p_scaler_id_zdc_wide   ,  "SCALER");
    long long unsigned int gl1p_scaler_zdc_narrow = p_gl1p->iValue(gl1p_scaler_id_zdc_narrow ,  "SCALER");

    // Get TrigLive, Trigraw, Trigscaled
    // These are identical to DST values, but we do not use for now.
    unsigned int trigraw    = p_gl1->iValue(0,RAWTRIG);
    unsigned int trigscaled = p_gl1->iValue(0,SCALEDTRIG);
    unsigned int triglive   = p_gl1->iValue(0,LIVETRIG);

    // Get BBC Z Vertex - found in online monitoring
    // This does not give the same value of the DST, we can't use it.
    float bbc_zvtx = p_bbcll1->iValue(0,"VERTEX");

    // Get ZDC Z Vertex

    /* 
    std::cout << "Event: " << event_number 
	    << ", raw: " << trigraw 
	    << ", live: " << triglive 
	    << ", scaled: " << trigscaled 
	    << ", bbcz: " << bbc_zvtx << std::endl;
	    */
    
    /* Debugging Output 
    std::cout << "EVENT NUMBER: " << event_number << std::endl;
    std::cout << "ATP#: " << atp_number << std::endl; // ATP done
    std::cout << "TIME: " << time << std::endl;
    std::cout << "BUNCH NUMBER: " << gl1_crossing_id << std::endl;
    std::cout << "GL1P_SCALER_BBC        : " << gl1p_scaler_bbc        << std::endl;
    std::cout << "GL1P_SCALER_CLOCK      : " << gl1p_scaler_clock      << std::endl;
    std::cout << "GL1P_SCALER_ZDC_WIDE   : " << gl1p_scaler_zdc_wide   << std::endl;
    std::cout << "GL1P_SCALER_ZDC_NARROW : " << gl1p_scaler_zdc_narrow << std::endl;
    */

    // Machine Readable Output
    out_file << event_number 
      << " " << atp_number 
      << " " << time 
      << " " << gl1_crossing_id 
      << " " << gl1p_scaler_bbc 
      << " " << gl1p_scaler_clock 
      << " " << gl1p_scaler_zdc_wide 
      << " " << gl1p_scaler_zdc_narrow 
      << std::endl;
  }
  return 0;
}
