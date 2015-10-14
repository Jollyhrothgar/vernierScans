#ifndef __HOURGLASS_SIMULATION__
#define __HOURGLASS_SIMULATION__

#include <string>
#include <vector>
#include "TH1F.h"
#include "TObject.h"

class HourglassSimulation {
  public:
   std::string this_name; // Stores the name of this class, plus it address in memory.
   HourglassSimulation();
   ~HourglassSimulation();
   int Init();
   int Run();  /** Simulates Z-vertex profile for ZDC, given configuration loaded. */
   int SaveFigures( const std::string& figure_output_dir); /** saves any figures loaded into save_restitry_ to a root file and pdf */
  private:
   int Factorial( int n ); /** compute factorial of integer n recursively */
   double SmearZVertex( double rand_prob_res, double orig_z ); /** given random probabiltiy and ZDC resolution, smear z-vertex */

   int SaveFigures();
   TH1F* zvertex;
   std::vector<TObject*> save_registry_;

   /** CONFIGURATION (These variables should not change after initialization...) )*/
   std::string run_number_;
  double xoff; 
  double yoff;  
  int count_norm;
  double rate; /**  RATE OF COLLISION AT EACH CROSSING (typical Run 6 'rate' ~ 0.129 */
  int max_coll; /** maximum collision number in a cross. = 5 */

  /** variables for z-t ~Gaussian distribution of intial beam profile */
  double n_bunch; /**  number of good bunches */
  double freq; /** freq. of bunch crossing */
  double scale;
  double Np; /**  run 12, proton density in a bunch is ~ 10^9 */
  double Nm; /**  run 12 */
  double sigma_xstar; /** cm */
  double sigma_ystar; /**  cm */
  double sigma_zl; /**  52.7 width of left Gaussian //cm */
  double sigma_zc; /**  41.48 width of center Gaussian //cm */
  double sigma_zr;  /**  width of right Gaussian  //cm */
  double mu_zl; /**  cm wp what is this parameter */
  double mu_zr; /**  cm wp what is this parameter */
  double  z_vtx_off; /**  offset from data for zdc zvtx, cm */
  double beta_star; /**  cm */
  double angle; /**  crossing angle(-0.2<->0.2) mrad  */
  
  /** Defining spatial corrdinates, to be used in arrays */
  double z_low; /** total length of beam ~ 600 cm */
  double z_high;
  double z_range;
  double x_low; /** total spread(20 sigma) of beam ~ 6 mm */
  double x_high;
  double x_range;
  double y_low; /** total spread(20 sigma) of beam ~ 6 mm */
  double y_high;
  double y_range;
  double t_low; /** total interaction time ~ 20 nano-sec */
  double t_high;
  double t_range;
  double vel; /**  speed of light (cm/s) */
  int N_bin_t;
  int N_bin_z;
  int N_bin_x;
  int N_bin_y;
  double binsizeG;
  double binsizeX;
  double binsizeY;
  double binsizeT;
  /** End spatial discreetization */
};

#endif
