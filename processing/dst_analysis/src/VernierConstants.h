#ifndef _VERNIER_CONSTANTS_H_
#define _VERNIER_CONSTANTS_H_

#include <string>

const std::string base_dir     = "/direct/phenix+spin2/beaumim/vernierScans/";
const std::string analysis_dir = base_dir+"vernier_analysis/";
const std::string figures_dir  = analysis_dir+"figs/";
const std::string macros_dir   = analysis_dir+"macros/";
const std::string data_dir     = base_dir+"data/";
const std::string trees_dir    = data_dir+"run_12/";
const std::string lists_dir    = analysis_dir+"lists/";
const std::string vernierTreeName = "reduced_vernier_DST";

const int NUMBER_OF_BUNCHES = 120;

#endif
