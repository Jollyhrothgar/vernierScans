#include <iostream>
#include <fstream>
#include <map>
#include <sstream>

struct merged_data{
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

int main(int argc, char** argv) {
	if(argc != 4){
    std::cout << " Usage is: " 
        << argv[0] << " <dst_file> <prdf_file> <output_file> " << std::endl;
    return 1;
  } 
  std::string dst_file = argv[1];
  std::string prdf_file = argv[2];
  std::string output_file = argv[3];
  std::cout << " Merging " << dst_file << " with " << prdf_file << std::endl;
  std::cout << " Sending output to : " << output_file << std::endl;

  std::map<long int, merged_data> final;
  std::map<long int, merged_data> dst;

  std::ifstream in_dst(dst_file.c_str());
  std::ifstream in_prdf(prdf_file.c_str());

  if(!in_dst) {
    std::cout << " Could not open " << dst_file << std::endl;
    return 1;
  }
  if(!in_prdf) {
    std::cout << " Could not open " << prdf_file << std::endl;
    return 1;
  }

  std::string line;
  std::cout << "  Reading DST" << std::endl;
  while(getline(in_dst,line)){
    merged_data d;
    d.Reset();
    long int event_number = 0;
    std::stringstream ss;
    ss << line;
    ss >> event_number 
       >> d.trigraw
       >> d.triglive
       >> d.trigscaled
       >> d.bbc_z 
       >> d.zdc_z;
    dst[event_number] = d;
  }
  std::cout << "  Loaded " << dst.size() << " entries from the dst file. " << std::endl;
  
  std::cout << "  Reading and merging with PRDF" << std::endl;
  while(getline(in_prdf,line)){
    merged_data d;
    d.Reset();
    long int event_number = 0;
    std::stringstream ss;
    ss << line;
    ss >> event_number 
      >> d.atp_number 
      >> d.timestamp 
      >> d.bunch 
      >> d.gl1p_bbc 
      >> d.gl1p_clock 
      >> d.gl1p_zdc_wide 
      >> d.gl1p_zdc_narrow ;
    auto i = dst.find(event_number);
    if(i != final.end()) { // Merge in data from DST
     d.trigraw     = (i->second).trigraw   ;
     d.triglive    = (i->second).triglive  ;
     d.trigscaled  = (i->second).trigscaled;
     d.bbc_z       = (i->second).bbc_z     ;
     d.zdc_z       = (i->second).zdc_z     ;
    }
    final[event_number] = d;
  }
  std::cout << " Merge complete. Final data set: " << final.size() << " entries from the prdf file. " << std::endl;
  std::ofstream out_file(output_file.c_str());
  for(auto i = final.begin(); i != final.end(); ++i ) {
    long int event_number = i->first;
    merged_data d = i->second;
    out_file << event_number 
      << " " <<  d.atp_number 
      << " " <<  d.timestamp 
      << " " <<  d.bunch 
      << " " <<  d.gl1p_bbc 
      << " " <<  d.gl1p_clock 
      << " " <<  d.gl1p_zdc_wide 
      << " " <<  d.gl1p_zdc_narrow 
      << " " <<  d.trigraw    
      << " " <<  d.triglive   
      << " " <<  d.trigscaled 
      << " " <<  d.bbc_z      
      << " " <<  d.zdc_z      
      << " " <<  std::endl;
  }
  return 0;
}
