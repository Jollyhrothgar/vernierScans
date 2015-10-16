#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

int main(int argc, char** argv) {
  if(argc < 2) {
    std::cout << "You need to use like so: " << std::endl
        << argv[0] << " file_to_count" << std::endl;
    return 1;
  }
  std::string in_file_name = argv[1];
  std::ifstream in_file(in_file_name.c_str());
  std::string line = "";
  long long unsigned int scaler_0_sum = 0;
  long long unsigned int scaler_1_sum = 0;
  long long unsigned int scaler_2_sum = 0;
  long long unsigned int scaler_3_sum = 0;

  if(in_file) {
    while(getline(in_file,line)){
      if(line[0] == '#') continue;
      std::stringstream ss;
      ss << line;
      long unsigned int evt_seq, atp_num, time, scaler_0, scaler_1, scaler_2, scaler_3, bunch_num;
      ss >> evt_seq >> atp_num >> time >> bunch_num >> scaler_0 >> scaler_1 >> scaler_2 >> scaler_3;
      scaler_0_sum += scaler_0;
      scaler_1_sum += scaler_1;
      scaler_2_sum += scaler_2;
      scaler_3_sum += scaler_3;
    }
  } else {
    std::cout << "Coule not open file: " << in_file_name << std::endl;
    return 1;
  }
  std::cout << "GL1P SCALER 0 SUM: " << scaler_0_sum << std::endl;
  std::cout << "GL1P SCALER 1 SUM: " << scaler_1_sum << std::endl;
  std::cout << "GL1P SCALER 2 SUM: " << scaler_2_sum << std::endl;
  std::cout << "GL1P SCALER 3 SUM: " << scaler_3_sum << std::endl;
  return 0;
}
