# Luminosity Calculation

Note: to view math properly, the extension 'tex everywhere'

This project template was generated with the module:
https://github.com/Jollyhrothgar/root_project_generator

The module consists of two classes, Luminosity, and CrossSection. These classes
are in charge of calculating Luminosity and Cross Section. It assumes that we
have already produced text files containing the relevant data. The data is
generated with the modules in
https://github.com/Jollyhrothgar/vernierScans/processes. 

# Equations

\begin{equation}
\mathcal{L} = f_{bunch} K_{\beta^*}K_{\theta_{xing}}\sum_{i=0}^{N_{bunches}}
{{N_{b_i}N_{y_i}}\over{2\pi\sigma_x\sigma_y}}
\end{equation}

\begin{equation}
\sigma_{BBC} = {{R_{BBC}}\over{\epsilon_{BBC}\mathcal{L}}}
\end{equation}

# Description of the Required data

## Constants

This analysis assumes a fixed bunch crossing frequency of 9.36 MHz.

## Beam Width Data

All bunches are assumed to have the same beam width in the transverse
directions. Although single bunches can be used to calculate the beam width,
some bunches have low statistics, and do not yield a suitable fit. Therefore,
we integrate the data for all the bunches. For a particular PHENIX run (say
359711), need a file, produced by:

vernierScans/processing/dst_analysis/macros/Run_BeamWidth.C

Which contains:

* Maximum overlap BBC Rate
* Horizontal Beam Width
* Vertical Beam Width

## Corrective Factors

Corrective factors are obtained by finding the optimum values for 
