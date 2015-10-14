#ifndef __BEAMSEPARATIONDATA__
#define __BEAMSEPARATIONDATA__
struct BeamSeparationData{
  double x; /** horizontal separation */
  double y; /** vertical separation */ 
  double x_err; /* uncertainty in x */
  double y_err; /* uncertainty in y */

  void Reset(){
    x = 0.;
    y = 0.;
    x_err = 0.;
    y_err = 0.;
  }
};
#endif
