#include <iostream>
#include <pmonitor.h> // linked to pmonitor libraries look for the root functions, prun, pinit, etc.
#include <EventTypes.h> // Added because it contains the definition for the Event ID in plain-english. Check out $ONLINE_MAIN/includes/

#include <fstream>
#include "clock.h"


int init_done = 0;

//TH1F *h1; 
//TH1F *h2; 

using namespace std;
static std::ofstream out_file;
// Call before running, it will do the right thing.
int set_output_file(const std::string& out_file_name) {
  out_file.open(out_file_name.c_str());
  out_file << "# epoch_time "  
      << " " << "clock_scaler_raw"
      << " " << "bbc_w_scaler_raw"
      << " " << "bbc_n_scaler_raw"
      << " " << "zdc_w_scaler_raw"
      << " " << "clock_scaler_live"
      << " " << "bbc_w_scaler_live"
      << " " << "bbc_n_scaler_live"
      << " " << "zdc_w_scaler_live"
      << endl;
  return 0;
}
// Call to gracefully exit
int finish() {
  out_file.close();
  return 0;
}

int pinit()
{

  if (init_done) return 1;
  init_done = 1;

  //  h1 = new TH1F("H1", "Channel 0", 101,-49.5,49.5);

  return 0;
}

// This is called by "prun", "pstart", etc.
int process_event (Event * e)
{
  // Run 12 Bits
  const int CLOCK_BIT = 11; // this is "Clock"
  const int BBC_NARROW = 0; // this is "BBCLL1(>0 tubes)" 
  const int BBC_WIDE = 1;   // this is "BBCLL1(>0 tubes) novertex"
  const int ZDC_WIDE = 2;   // this is "ZDCLL1wide"
  if(e->getEvtType() == SCALEREVENT) // SCALEREVENT is event type 14
  {
    //e->identify(); // prints out general info about the packet.
    // Static insures this variable's value is stored between calls to 
    // this function
    static unsigned long long clock_scaler_raw_old = 0;
    static unsigned long long bbc_w_scaler_raw_old = 0;
    static unsigned long long bbc_n_scaler_raw_old = 0;
    static unsigned long long zdc_w_scaler_raw_old = 0;
    static unsigned long long clock_scaler_live_old = 0;
    static unsigned long long bbc_w_scaler_live_old = 0;
    static unsigned long long bbc_n_scaler_live_old = 0;
    static unsigned long long zdc_w_scaler_live_old = 0;
    unsigned long long clock_scaler_raw_current;
    unsigned long long bbc_w_scaler_raw_current;
    unsigned long long bbc_n_scaler_raw_current;
    unsigned long long zdc_w_scaler_raw_current;
    unsigned long long clock_scaler_live_current;
    unsigned long long bbc_w_scaler_live_current;
    unsigned long long bbc_n_scaler_live_current;
    unsigned long long zdc_w_scaler_live_current;
    static unsigned long long clock_scaler_raw_offset = 0;
    static unsigned long long bbc_w_scaler_raw_offset = 0;
    static unsigned long long bbc_n_scaler_raw_offset = 0;
    static unsigned long long zdc_w_scaler_raw_offset = 0;
    static unsigned long long clock_scaler_live_offset = 0;
    static unsigned long long bbc_w_scaler_live_offset = 0;
    static unsigned long long bbc_n_scaler_live_offset = 0;
    static unsigned long long zdc_w_scaler_live_offset = 0;
    
    //cout << sizeof(zdc_w_scaler_raw_current) << endl;
    Packet *p = e->getPacket(900); // scaler packet id
    if (p)
    {
      // h1->Fill (p->iValue(0));
      // dlist to get a list of packets from the command line
      // the identifier for the "RAWSCALER" type is the scaler bit (i.e. 0-31) 
      //cout << p->iValue(17,"RAWSCALERS") << " " << e->getTime() << endl;
      // These values may have overflowed (clock, bbcwide)
      clock_scaler_raw_current = (unsigned int)p->iValue(CLOCK_BIT,"RAWSCALERS") ;
      bbc_w_scaler_raw_current = (unsigned int)p->iValue(BBC_WIDE,"RAWSCALERS")    ; 
      bbc_n_scaler_raw_current = (unsigned int)p->iValue(BBC_NARROW,"RAWSCALERS")   ;
      zdc_w_scaler_raw_current = (unsigned int)p->iValue(ZDC_WIDE,"RAWSCALERS")     ;
      clock_scaler_live_current = (unsigned int)p->iValue(CLOCK_BIT,"LIVESCALERS") ;
      bbc_w_scaler_live_current = (unsigned int)p->iValue(BBC_WIDE,"LIVESCALERS")    ; 
      bbc_n_scaler_live_current = (unsigned int)p->iValue(BBC_NARROW,"LIVESCALERS")   ;
      zdc_w_scaler_live_current = (unsigned int)p->iValue(ZDC_WIDE,"LIVESCALERS")     ;

      // We check if the most significant bit has changed here, in order to deterine if
      // an overflow has occurred. We keep track of the the extra bit, and add it to our total.
      if(clock_scaler_raw_current < clock_scaler_raw_old) {
        clock_scaler_raw_offset += 0x100000000ULL;
      }
      if(bbc_w_scaler_raw_current < bbc_w_scaler_raw_old) {
        bbc_w_scaler_raw_offset += 0x100000000ULL;
      }
      if(bbc_n_scaler_raw_current < bbc_n_scaler_raw_old) {
        bbc_n_scaler_raw_offset += 0x100000000ULL;
      }
      if(zdc_w_scaler_raw_current < zdc_w_scaler_raw_old) {
        zdc_w_scaler_raw_offset += 0x100000000ULL;
      }
      if(clock_scaler_live_current < clock_scaler_live_old) {
        clock_scaler_live_offset += 0x100000000ULL;
      }
      if(bbc_w_scaler_live_current < bbc_w_scaler_live_old) {
        bbc_w_scaler_live_offset += 0x100000000ULL;
      }
      if(bbc_n_scaler_live_current < bbc_n_scaler_live_old) {
        bbc_n_scaler_live_offset += 0x100000000ULL;
      }
      if(zdc_w_scaler_live_current < zdc_w_scaler_live_old) {
        zdc_w_scaler_live_offset += 0x100000000ULL;
      }
      // Check that the overflow correction is working.
      // cout << e->getTime() 
      //     << " " << (unsigned int)p->iValue(CLOCK_BIT,"RAWSCALERS")  << ":" << clock_scaler_raw_current + clock_scaler_raw_offset
      //     << " " << (unsigned int)p->iValue(BBC_WIDE,"RAWSCALERS")   << ":" << bbc_w_scaler_raw_current + bbc_w_scaler_raw_offset
      //     << " " << (unsigned int)p->iValue(BBC_NARROW,"RAWSCALERS") << ":" << bbc_n_scaler_raw_current + bbc_n_scaler_raw_offset
      //     << " " << (unsigned int)p->iValue(ZDC_WIDE,"RAWSCALERS")   << ":" << zdc_w_scaler_raw_current + zdc_w_scaler_raw_offset
      //     << endl;
      out_file << e->getTime() 
          << " " << clock_scaler_raw_current + clock_scaler_raw_offset
          << " " << bbc_w_scaler_raw_current + bbc_w_scaler_raw_offset
          << " " << bbc_n_scaler_raw_current + bbc_n_scaler_raw_offset
          << " " << zdc_w_scaler_raw_current + zdc_w_scaler_raw_offset
          << " " << clock_scaler_live_current + clock_scaler_live_offset
          << " " << bbc_w_scaler_live_current + bbc_w_scaler_live_offset
          << " " << bbc_n_scaler_live_current + bbc_n_scaler_live_offset
          << " " << zdc_w_scaler_live_current + zdc_w_scaler_live_offset
          << endl;

      // Raw Scalers
      clock_scaler_raw_old = clock_scaler_raw_current;
      bbc_w_scaler_raw_old = bbc_w_scaler_raw_current;
      bbc_n_scaler_raw_old = bbc_n_scaler_raw_current;
      zdc_w_scaler_raw_old = zdc_w_scaler_raw_current;

      // Live Scalers
      clock_scaler_live_old = clock_scaler_live_current;
      bbc_w_scaler_live_old = bbc_w_scaler_live_current;
      bbc_n_scaler_live_old = bbc_n_scaler_live_current;
      zdc_w_scaler_live_old = zdc_w_scaler_live_current;

      delete p;
    }
  }
  return 0;
}

