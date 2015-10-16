#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <fstream>

struct run_data {
  void Reset() {
    gl1p_1 = 0;
    gl1p_2 = 0;
    gl1p_3 = 0;
    gl1p_4 = 0;
  }
  long int gl1p_1;
  long int gl1p_2;
  long int gl1p_3;
  long int gl1p_4;
};

int main(int argc, const char ** argv) {
  if( argc != 4 ) {
    std::cout << "Usage is: " << argv[0] 
        << " input_file output_prefix output_directory" 
        << std::endl;
    return 1;
  }

  std::map < long int, run_data > run_integrated;
  std::map < int , std::map<long int, run_data > > bunch_data;
  for( int i = 0; i < 120; i++ ) {
    bunch_data[i];
  }


  std::string input_file = argv[1];
  std::string output_prefix = argv[2];
  std::string output_directory = argv[3];
  
  std::cout << "Reading " << input_file << std::endl
      << "   labeling output with " << output_prefix << std::endl
      << "   sending output to " << output_directory
      << std::endl;
  std::ifstream in_file(input_file.c_str());
  std::string line = "";
  if(in_file) {
    while(getline(in_file,line)) {
      std::stringstream ss;
      long int time_stamp;
      int atp_num;
      long int evt_num;
      int bunch_number;
      long int gl1p_1;
      long int gl1p_2;
      long int gl1p_3;
      long int gl1p_4;
      ss << line;
      ss 
        >> evt_num 
        >> atp_num 
        >> time_stamp 
        >> bunch_number 
        >> gl1p_1
        >> gl1p_2
        >> gl1p_3 
        >> gl1p_4;
      if ( run_integrated.find(time_stamp) == run_integrated.end() ) {
        // not found, need to add
        run_integrated[time_stamp];
        run_integrated[time_stamp].Reset();
        run_integrated[time_stamp].gl1p_1 = gl1p_1;
        run_integrated[time_stamp].gl1p_2 = gl1p_2;
        run_integrated[time_stamp].gl1p_3 = gl1p_3;
        run_integrated[time_stamp].gl1p_4 = gl1p_4;
        bunch_data[bunch_number][time_stamp];
        bunch_data[bunch_number][time_stamp].Reset();
      } else {
        // found, index count
        run_integrated[time_stamp].gl1p_1 += gl1p_1;
        run_integrated[time_stamp].gl1p_2 += gl1p_2;
        run_integrated[time_stamp].gl1p_3 += gl1p_3;
        run_integrated[time_stamp].gl1p_4 += gl1p_4;
        bunch_data[bunch_number][time_stamp].gl1p_1 += gl1p_1;
        bunch_data[bunch_number][time_stamp].gl1p_2 += gl1p_2;
        bunch_data[bunch_number][time_stamp].gl1p_3 += gl1p_3;
        bunch_data[bunch_number][time_stamp].gl1p_4 += gl1p_4;
      }
    }
  } else {
    std::cout << "could not open file: " << input_file << std::endl;
    return 1;
  }

  std::string out_file_name = output_directory + "/" + output_prefix + "_bunch_integrated.txt";
  std::ofstream out_file(out_file_name.c_str());
  for(auto i = run_integrated.begin(); i != run_integrated.end(); ++i) {
    out_file << i->first
      << " " << i->second.gl1p_1
      << " " << i->second.gl1p_2
      << " " << i->second.gl1p_3
      << " " << i->second.gl1p_4
      << std::endl;
  }
  out_file.close();
  for(int bunch_i = 0; bunch_i < 120; bunch_i++) {
    std::stringstream bunch_file_name;
    bunch_file_name << output_directory << "/" << output_prefix + "_bunch_" << bunch_i << ".txt";
    std::ofstream bunch_file(bunch_file_name.str().c_str());
    for( auto i = bunch_data[bunch_i].begin(); i != bunch_data[bunch_i].end(); ++i) {
      bunch_file << i->first
        << " " << i->second.gl1p_1
        << " " << i->second.gl1p_2
        << " " << i->second.gl1p_3
        << " " << i->second.gl1p_4
        << std::endl;
    }
  }
  return 0;
}
