#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
/** Instructions
 * This is a command line utility. It is desined to work with 
 * raw prdf data which is dumped to a text file. Each line of data 
 * dumped from a prdf to a text file must contain a minimum of a
 * unix time stamp, and an atp number. This code takes a list of 
 * dumped prdf segments as an argument, corrects the time offset
 * on the unix time stamp for each atp, and dumps the corrected data
 * out to files in a different directory.
 *
 * This utility does not have any functionality in terms of guessing
 * or being told what order your data is in, so you have to explitly
 * hard code it in, or, write a better utility, =).
 *
 * This ulility assumes that the events in each dumped prdf are 
 * in the same order as the events that were recorded in the prdf, 
 * and exploits this to assume rough time ordering. If you mess with
 * the ordering of your data,  this utiltiy will break. If you 
 * merge your segments before running over your data, this 
 * utility will break.
 *
 * Basic work flow:
 * 1) Dump each PRDF segment to its own text file
 * 2) create a file containing the full path to each segment + the
 *    segment's name, with one /<path>/seg_name.txt per line.
 * 3) Compile this utility, such as ":~> g++ correct_time.C -o correct_time.exe"
 * 4) Run the utility: ":~> correct_time.exe file_list.txt merged_corrected_data.txt"
 * Author: Mike Beaumier
 * Contac: michael.beaumier@gmail.com
 */
int main (int argc, char** argv){

	if(argc != 3) {
		std::cout << "example usage is " << argv[0] << "run_list.txt merged_data.txt"  << std::endl;
		return -1;
	}

	// setup each ATP and initialize
	int atptime[64];
	int atpcount[64];
	int atpcorr[64];

	for ( int i = 0; i < 64; i++){
		atptime[i] = 0;
		atpcount[i] = 0;
	}
	std::string atp_files = argv[1];
	std::cout << "  Loading PRDF dump list: " << atp_files << std::endl;
	std::ifstream in_list(atp_files.c_str());
	std::vector<std::string> file_list;
	std::string file = "";
	while(getline(in_list,file)) {
		std::cout << "    Got file name: " << atp_files << std::endl;
		file_list.push_back(file);
	}
	in_list.close();

	std::vector<std::string>::iterator file_i;
	for(file_i = file_list.begin(); file_i != file_list.end(); ++file_i) {
		int count = 0;
		std::cout << "Getting time offset for file: " << *file_i << std::endl;
		std::ifstream in_file((*file_i).c_str());

		if(in_file) {
			std::string line = "";
			while(getline(in_file,line) && count++ < 20000) {
				std::stringstream ss;
				long int time_stamp;
				int atp_num;
				long int evt_num;
				ss.str(line);
				ss >> evt_num >> atp_num >> time_stamp; 
				// only do a partial extraction here, we don't need the rest yet
				// Here, we just find the first instance of the ATP being used.
				if ( atp_num < 64 && atpcount[atp_num] == 0 ) {
					atptime[ atp_num ] +=  time_stamp ;
					atpcount[ atp_num ]++;
				}
			}
			std::cout << "  done" << std::endl;
		}
		else {
			std::cout << "can't open " << *file_i << std::endl;    
		}
	}
	std::cout << "ATP OFFSET: " << std::endl;
	std::cout << "      " 
		<< std::setw(3) <<  "atp" 
		<< std::setw(12) << "time_stamp" 
		<< std::setw(12) << "offset" 
		<< std::setw(12) << "corr" << std::endl;
	for (int atp_nr = 0; atp_nr < 64; atp_nr++) {
		if  (atpcount[atp_nr] ) {
			atpcorr[atp_nr] = atptime[atp_nr]-atptime[0]; // synchronize to atp 0
			std::cout 
				<< "      " << std::setw(3) <<  atp_nr 
				<< std::setw(12) << atptime[atp_nr] 
				<< std::setw(12) << atpcorr[atp_nr] 
				<< std::setw(12) << atptime[atp_nr] - atpcorr[atp_nr] 
				<<  std::endl;
		}
	}
	// Now, we may correct the time stamp for all ATPs.
	std::string out_file_name = argv[2];
	out_file_name = "./" + out_file_name;
	std::cout << "Sending corrected data to: " << out_file_name << std::endl;
	std::ofstream out_file(out_file_name.c_str());
	for(file_i = file_list.begin(); file_i != file_list.end(); ++file_i) {
		std::cout << "Correcting : " << *file_i << std::endl;
		int count = 0;
		std::ifstream in_file((*file_i).c_str());

		if(in_file) {
			std::string line = "";
			while(getline(in_file,line)) {
				std::stringstream ss;
				// This may need to be changed!!
				// 8074318 7 1329589837 31 2 71 0 0
				long int time_stamp;
				int atp_num;
				long int evt_num;
				int bunch_number;
				long int gl1p_1;
				long int gl1p_2;
				long int gl1p_3;
				long int gl1p_4;

				ss.str(line);
				ss 
					>> evt_num 
					>> atp_num 
					>> time_stamp 
					>> bunch_number 
					>> gl1p_1
					>> gl1p_2
					>> gl1p_3 
					>> gl1p_4;
				//correct:
				if(atp_num >= 0 && atp_num < 64) {
					time_stamp = time_stamp - atpcorr[atp_num];
					// dump
					out_file  
						<< evt_num 
						<< " " << atp_num 
						<< " " << time_stamp 
						<< " " << bunch_number 
						<< " " << gl1p_1
						<< " " << gl1p_2
						<< " " << gl1p_3 
						<< " " << gl1p_4
						<< std::endl;
				}
			}
		} else {
			std::cout << "can't open " << *file_i << std::endl;    
		}
		std::cout << " Done." << std::endl;
	}
	out_file.close();
	std::cout << "Everything is finished!" << std::endl;
	return 0;
}
