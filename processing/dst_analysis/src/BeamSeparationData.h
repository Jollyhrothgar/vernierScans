#ifndef __BEAMSEPARATIONDATA__
#define __BEAMSEPARATIONDATA__
struct BeamSeparationData{
  double x; /** horizontal separation */
  double y; /** vertical separation */ 
  double x_err; /* associated uncertainty in x */
  double y_err; /* associated uncertainty in y */

  double x_avg; /* horizontal separation calculated from average beam pos */
  double y_avg; /* vertical separation calculated from avg beam pos */
  double x_avg_err; /* associated uncertainty of x_avg */
  double y_avg_err; /* associated uncertainty of y_avg */


  void Reset(){
    x = 0.;
    y = 0.;
    x_err = 0.;
    y_err = 0.;
  }
};
#endif
