#ifndef __BEAM_WIDTH_DATA__
#define __BEAM_WIDTH_DATA__

#include <cmath>

struct BeamWidthData {
  BeamWidthData() { Reset();};
  double rate; /** scaler rate (BBC) */
  double rate_stat_err;
  double rate_sys_err;
  double x; /** Horizontal beam separation */
  double y; /** Vertical beam separation */
  double x_planned; /** Planned horizontal separation */
  double y_planned; /** Planned vertical separation */
  double x_stat_err;
  double y_stat_err;
  double x_sys_err;
  double y_sys_err;
  void Reset() {
    rate = 0.;
    rate_stat_err = 0.;
    rate_sys_err = 0.;
    x_planned = 0.;
    y_planned = 0.;
    x = 0.;
    y = 0.;
    x_stat_err = 0.;
    y_stat_err = 0.;
    x_sys_err = 0.;
    y_sys_err = 0.;
  };
  double GetTotalErrX(){
    return pow(pow(x_sys_err,2.0)+pow(x_stat_err,2.0),0.5);
  }
  double GetTotalErrY(){
    return pow(pow(x_sys_err,2.0)+pow(x_stat_err,2.0),0.5);
  }
  double GetTotalErrRate() {
    return pow(pow(rate_stat_err,2.0)+pow(rate_sys_err,2.0),0.5);
  }
};
#endif
