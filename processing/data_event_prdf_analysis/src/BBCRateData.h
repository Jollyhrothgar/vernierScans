#ifndef __BBCRATEDATA__
#define __BBCRATEDATA__
/** This is meant to be a container, though we could add
 * various member functions to calculate rate or rate_err
 */
struct BBCRateData {
  double rate;         /** BBC Rate, want from trigger: BBCLL1(>0 tubes) */
  double rate_err;     /** statistical error of the rate */                
  long int bbc_gl1p;   /** # of bbc gl1p counts used to calculate rate */  
  long int clock_gl1p; /** # of clock gl1p counts used to calculate rate */
  void Reset() {
    rate = 0.;
    rate_err = 0.;  
    bbc_gl1p = 0;   
    clock_gl1p = 0; 
  }
};
#endif
