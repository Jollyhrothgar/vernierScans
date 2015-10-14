#ifndef _VERNIER_TREE_VARIABLES_H_
#define _VERNIER_TREE_VARIABLES_H_

#include "TTree.h"
#include <map>
#include <string>
#include <iostream>

struct VernierTreeVariables
{
  VernierTreeVariables();
  int   cross_id     ;
  int   evtnumber    ;
  int   trigraw      ;
  int   triglive     ;
  int   trigscaled   ;
  double bbc_z        ;
  double zdcll1_z     ;
  int   BbcNorthTubes;
  int   BbcSouthTubes;
  
  std::map<std::string, int> choice;

  int Reset();

  int LinkTree(TTree*, const std::string&);
};
#endif
