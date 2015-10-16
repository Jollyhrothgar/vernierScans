#ifndef __RUN_SCALERS_H__
#define __RUN_SCALERS_H__

#include <Event.h>

int process_event (Event *e);

int set_output_file(const std::string& out_file);

int finish(); 

#endif /* __RUN_SCALERS_H__ */
