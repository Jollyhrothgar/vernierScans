This is the complete vernier analysis software repository written by Mike Beaumier. The
'processing' portion of the repository is self-contained, and uses data sorces including
root files (containing root trees and histograms) and text files. The "reduce_data"
portion of the repository may use condor job management, proprietary PHENIX libraries, or
other tools which may only be accessed within the RCF computing environment.

Both sub-respositories depend on the "FileManagement.h" header, in this directory, which
allows for easier job management - various submodules use compiled libraries which are
called in root macros. These root macros often take data, organized into arrays of various
file-handles in FileManagement.h, such that you may run an analysis on a single Vernier
scan by simply choosing the appropriate array index, shared by all the containers in
FileManagement.h.

reduce_data - chunks down the raw data sources, applying corrections when neccessary.
processing - runs the analysis over the output data of reduce_data

both depend on data stored in a local folder here, "data" with the directory structure:

  vernierScans/data/run_12/DSTs               -> Raw DSTs from Fun4All vernier production 
  vernierScans/data/run_12/prdf_data          -> Corrected, dumped data-events from PRDFF
  vernierScans/data/run_12/scaler_data        -> Scaler events from ppg PRDFF data
  vernierScans/data/run_12/prdf_dst_merged    -> Merged output output of dst_data and
                                                 prdf_data 
  vernierScans/data/run_12/simulation_config  -> Configuration files used by
                                                 HourglassCorrection simulation 
  vernierScans/data/run_12/cad_data           -> BPM, WCM, DCCT data
  vernierScans/data/run_12/dst_data           -> Root trees produced by
                                                 verneirScans/reduce_data/dst_reduction
The overall workflow is kept modular, which is likely confusing to someone who is not
familiar with PHENIX software, or not familiar with the workflow of the vernier analysis.
I chose to structure my code modularly so that working on one part of the analysis will
not destroy all parts of the analysis. I'm not sure who I'm writing this readme for, since
I'm likely going to be the only one whoever uses this code. Practice I guess?
