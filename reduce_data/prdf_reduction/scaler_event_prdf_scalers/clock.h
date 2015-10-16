#ifndef __CLOCK_H__
#define __CLOCK_H__

#include <Event.h>
int process_event (Event *e);

int set_output_file(const std::string& out_file);

int finish(); 
#endif /* __CLOCK_H__ */
