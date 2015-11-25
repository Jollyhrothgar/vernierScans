#ifndef __HOURGLASS_SIMULATION__
#define __HOURGLASS_SIMULATION__

#include "HourglassConfiguration.h"
#include <string>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TObject.h"


class HourglassSimulation {
 public:
  // Stores the name of this class, plus it address in memory.
  std::string this_name; 
  
  HourglassSimulation();
  ~HourglassSimulation();
  
  // Skip the config file, and initalize everything from hardcoded parameters.
  // This is for emergency debugging, I wouldn't reccomend using it for anything
  // else.
  int ManualInit();

  // Use HourglassConfiguration to generate and load the default configuration
  // which is done manually in InitDefault. InitDefault is kept as a debugging
  // tool.
  int InitDefaultConfig();

  // Loads all configuration parameters from HourglassConfig-generated file, use
  // for running many simulations over smoothly varying parameters.
  int InitFromConfig(const std::string& config_file);

  // For Development/shape Tuning: override the default save file location and
  // name to run simulaitons interactively in order to tune simulation.
  int OverrideSaveFile(const std::string& file_name);
  std::string override_save_file_name_;
  bool override_save_file_;

  // Simulates Z-vertex profile for ZDC, given configuration loaded. 
  int Run(int model_opt);  

  // Shows the simulation configuration used. Some config parameters have
  // operations applied to them immediately before use.
  int ShowConfig();

  // Compare to a data distribution
  // The proper zdc distribution is extracted from the config file, the root
  // file which contains this distribution is passed as compare_file_name. The
  // simulated z-vertex is plotted on top of the real z-vertex distribution, and
  // the canvas is saved to the registry.
  int Compare(const std::string& compare_file_name);
  std::string zdc_compare_histo_name_;
  
  // saves any figures loaded into save_restitry_ to a root file and pdf. Used
  // to capture all graphical output and store it for inspection later, no
  // formatting is performed to make the plots pretty.
  int SaveFigures( const std::string& figure_output_dir); 
 private:
  // Helper function to pass the configuration parameters from an
  // HourglassConfiguration object to this object. 
  int InitConfig();

  // Initialized to the configuration file which was used to run this
  // simulation.
  HourglassConfiguration config_;

  // Initialize space-time related coordinates
  int InitSpacetime();

  // Performance Tracking Variables
  long long unsigned int how_many_things = 0; 
  std::map<int, long long unsigned int > time_tracker;

  // Spatial Coordinates
  std::vector<double> z_position_; 
  std::vector<double> x_position_; 
  std::vector<double> y_position_; 
  std::vector<double> t_position_;

  // Probability Distributions, and Related Variables
  std::vector<double> poisson_dist_;
  std::vector< double > t_dist_;
  std::vector< double > z_dist_;
  std::vector< double > z_norm_;
  std::vector< std::vector<double> > gaussian_dist_;
  unsigned long long cross_count;
  unsigned long long event_limit_count;

  // Luminosity variables, functions
  //
  // Find the normalization for our density, assumes un-normalized x and y
  // distributions, and samples the normalized z-distriubtion, which is already
  // normalized. Note that rotating the profile does not change the overall
  // normalization.
  void NormalizeDensity();

  // Standard gaussian, sans normalization
  double GetGaussianDensity(double x, double sigma);

  // Looks up z value and returns the density based on the zprofile loaded.
  // lookup value is first argument. The second argument updates z_lookup with
  // the value which was found in the distribution, so that z_lookup and z can
  // be compared.
  double LookupZDensityBlue(double z,double& z_lookup);
  double LookupZDensityYell(double z,double& z_lookup);

  // total luminosity calculated in numeric integration
  double luminosity_tot_; 
  double luminosity_normalization_;
  double density_normalization_;
  double spacetime_volume_;
  double x_profile_norm_;
  double y_profile_norm_;
  double z_profile_norm_;
  
  // Initializes plots, must be called after InitSpacetime, which is called by
  // either of the other init methods. I put it at the top of Run()
  int InitPlots();

  // Final initialization, once I figure out what is initialized here, I can
  // rename to something more sensible.
  int InitProbabilityVariables();

  // compute factorial of integer n recursively 
  int Factorial( int n ); 

  // The default model for luminosity we use. Only one crossing angle is
  // observable in the resultant z-vertex profile, assuming that the horizontal
  // and vertical scans are done separately (i.e. no diagonal scans).  Detects
  // the scan displacement direction from the configuration provided
  //
  // This function is relatively speedy, but requires that in our luminosity
  // calculation, that we have defined the following distriubtions and
  // variables:
  // 
  // * count_norm -> the number of ZDC counts
  // * gaussian_dist_ -> 
  // * z_norm_ -> 
  // * z_dist_ -> 
  int CreateCumulativePoissonDistribution();
  void GenerateAmareshModel();
  void GenerateNewModel();

  // Uses the output of our luminosity model to generate a z-vertex profile.
  // This is accomplished by sampling the z-t distribution created through
  // partial integration of the luminosity model.
  int GenerateZVertexProfile();

  // Assumes rotation in CZ plane, but rotation proceeds the same if we look at
  // the XZ plane or YZ plane. C = X, Y. Angle is not neccessarily the same in
  // both planes. Updates c and z with rotated values.
  //
  // It seems that rotations in three dimensional space are non commutative.
  // Which means, it matters what order we apply rotations in.
  void RotateBlueCoordinates(double& c, double& z, double angle);
  void RotateYellCoordinates(double& c, double& z, double angle);

  // given random probabiltiy and ZDC resolution, smear z-vertex 
  double SmearZVertex( double rand_prob_res, double orig_z ); 

  // Apply Beta Squeeze to beam_Width, return squeezed
  // beam width
  double BetaSqueeze(double beam_width, double z);

 public: 
  // Loading z-vertex distributions from WCM data. Distributions loaded to
  // z_bunch_profile_(blue|yell)_
  // Assumes the distribution is named as follows:
  // std::string b_name = "bwcm_zprofile_"+run_number_+"_space";
  // std::string y_name = "ywcm_zprofile_"+run_number_+"_space";
  // You can add an argument to this member function to make it general, but as
  // this is my code, and I'm trying to get results without messing around with
  // flexibility, the names are hardcoded.
  int LoadZProfile(const std::string& blue_f_name, const std::string& yell_f_name);
  std::map<double,double> z_profile_blue;
  std::map<double,double> z_profile_yell;
 
 private: 

  // Normalization:
  // In any luminosity model, we assume integration over the convelution of two
  // normalized densities. When those densities are simple, we know the
  // normalization ahead of time, but if we transform the widths and coordinates
  // to be multivariate distributions to add more "realism" the normalization
  // might not maintain the same functional form. Therefore, we should calculate
  // the normalization for each density individually, and apply it to the
  // luminosity convelution. This normalization should be trivial to calculate
  // numerically
  int NormalizeDensities(); 

  // Simulated z-vertex distribution
  TH1F* zdc_zvertex_sim;
  TH1F* zdc_zvertex_dat;
  TH1F* zdc_zvertex_sim_norm;
  TH1F* zdc_zvertex_dat_norm;
  TH2F* zvtx_pdf;
  TH1F* profile_frame_;
  TGraph* z_bunch_profile_blue_;
  TGraph* z_bunch_profile_yell_;
  TCanvas* zvertex_comparison_canvas;
  TCanvas* simulation_config_canvas;
  TCanvas* config_and_vertex_compare;
  TCanvas* norm_config_and_vertex_compare;
  std::vector<TObject*> save_registry_;

 public:
  // goodness of fits
  double chi2_test;
  double squares_residual;
 private:

  // CONFIGURATION (These variables should not change after initialization...) )
  std::string run_number_;

  // Beam Scan offsets
  double xoff, yoff; 
  int count_norm;
  //  Multiple collisions parameter, this is the typical number of collisions
  //  per bunch crossing. This varies with beam overlap, beam energy, beam
  //  population.
  double rate; 
  
  // Multiple collisions parameter, maximum expected number of collisions per
  // crossing is 0.5. MAX_COLL is a constant which is defined to be max_coll +
  // 1, and exists for the purpose of array indexing properly.
  int max_coll, MAX_COLL; 

  //  number of good bunch crossings
  double n_bunch; 

  // frequency of bunch crossings 
  double freq; 

  // scale the z-profile bunch geometry. Why we need this is not certain.
  double scale;
  
  // Two simulated bunches, N_blue and N_yell correspond to the average bunch
  // seen in the blue or yellow beam. Bunch asymmetries may lead to z-profile
  // skewing, so we need to be consistant. Unitless parameter, order of
  // magnitude is 10^9.
  double N_blue,N_yell; 

  // Beam widths obtained from vernier scan analysis. Horizontal scan width
  // (sigma_x) and vertical beam width (sigma_y). Because the luminosity
  // calculation defines sigma_x* and sigma_y* as sigma_x and sigma_y with an
  // additional factor of sqrt(2), we create these variables too.  unit: cm
  double sigma_xstar, sigma_ystar, sigma_x, sigma_y;
  
  // ZPROFILE BUNCH GEOMETRY:
  // Z-Profile of each simulated bunch is modeled with three gaussians, a
  // left-hand side gaussian, a right-hand side gaussian, and a central
  // gaussian. The parameters are as follows: left-hand gaussian width:
  // sigma_zl, right-hand gaussian width: sima_zr, center gausian width:
  // sigma_zc, left gaussian offset: mu_zl, right gaussian offset: mu_zr, with
  // the center gaussian assumed to have an offset of 0. Units are all cm.
  //
  // variables with sc_ prefix have had the "scale" member variable applied as a
  // mulitplicative factor.
  double sigma_zl, sigma_zr, sigma_zc, mu_zl, mu_zr, sc_sigma_zl, sc_sigma_zr,
         sc_sigma_zc, sc_mu_zl, sc_mu_zr;

  // Offset between BBC average z-vertex and ZDC average z-vertex, taken at
  // points where beams are maximally overlapped. Unit is cm.
  double  z_vtx_off; 
  
  // Beam focusing parameter, beta-star. Unit is in cm
  double beta_star; 
  
  // X-Z crossing angle(-0.2<->0.2) mrad  
  double angle_xz; 

  // Y-Z crossing angle(-0.2<->0.2) mrad
  double angle_yz;
  
  // Discreet Space-Time Variables
  
  // Z-coordinate of space is set to be approximately the length of one beam
  // bunch (600 cm).  N_bin_z determines the granualarity of the space-time
  // coordinate..
  double z_low, z_high, z_range; 
  int N_bin_z;
  double binsizeZ;
  
  // x-coordinate of space is set to be about 20 times larger than the beam
  // width (which is nominally 0.6 cm). N_bin_x determines the granualarity of
  // the space-time coordinate.
  double x_low, x_high, x_range;
  int N_bin_x;
  double binsizeX;
  
  // y-coordinate of space is set to be about 20 times larger than the beam
  // width (which is nominally 0.6 cm). N_bin_y determines the granualarity of
  // the space-time coordinate.
  double y_low, y_high, y_range;
  int N_bin_y;
  double binsizeY;
  
  // time-coordinate is set to be representative of the total interaction time,
  // which is about 20 nanoseconds. N_bin_t variable determines the granualarity
  // of the space-time coordinate.
  double t_low, t_high, t_range;
  int N_bin_t;
  double binsizeT;
  
  //  speed of light (cm/s) 
  double vel; 

  // End spatial discreetization 
};

#endif
