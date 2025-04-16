Monte Carlo Simulations for Estimation of IVIM Parameter Errors

This group of MATLAB scripts was written to support results provided in the supplemental materials of the whitepaper (publication pending) produced from the 2024 ISMRM Workshop on Moving Forward with Intravoxel Incoherent Motion Modeling for Diffusion-Weighted MRI, An Attempt at Consensus.

The intended use of these scripts is to allow users to estimate bias and precision errors for IVIM parameters for a selection of diffusion weightings (b-values) and signal-to-noise conditions. A brief description of these methods follow:

As outlined in the whitepaper, two ranges of b-values are recommended for consistency among the community, the so called “minimal” and “abbreviated”. Scripts to perform Monte Carlo simulations for these methods are named: DWI_IVIM_whitepaper_sim_minimal_bval.m and DWI_IVIM_whitepaper_sim_abbreviated_bval.m, respectively. 

From on the whitepaper, IVIM signal decay data are simulated with equation 1 using average IVIM parameters from the literature (Table 3) and both “minimal” and “abbreviated” b-value sets (Table 4). Simulations use complex data where the real channel included simulated IVIM decays combined with Gaussian additive noise and the imaginary channel included only Gaussian noise. Real and imaginary signals are combined as a root sum of squares and fit using nonlinear least squares. Signal to noise (SNR) is defined as the ratio of the b=0 signal amplitude divided by the standard deviation of the noise. 

For fitting of abbreviated b-value data, a monoexponential model is fit (IVIM_monoSEG_bifit.m) to decay data above the IVIM threshold value (see Table 4 and accompanying text) for each respective organ context to obtain D. The perfusion fraction, f, is estimated by subtracting the y-intercept obtained from the monoexponential fit from the b=0 signal intensity. 

For the minimal b-value data, a segmented fit is performed (IVIM_SEG_bifit_3par_opt1.m) where a monoexponential model is first fit to the decay data above the IVIM threshold value to obtain D and f as described in the minimal b-value fitting scheme using nonlinear least squares and the same fitting constraints. Next, the biexponential model is fit to all b-value data constraining the parameters f, 1-f, and D that were obtained in the first step while estimating D*. 

All fits are performed using nonlinear least squares with the following bounds on parameter estimates (see equations below): a=[0.1, 3]; D=[0.1e-3, 5e-3]; and D*=[0.5e-2, 2]

Bias and dispersion errors for each IVIM parameter are estimated for each SNR and organ context based on 1,000 independent noise realizations. Bias error is defined as the relative error 100*(estimate – true)/true and dispersion error is defined as the 100*coefficient of variation (CV = standard deviation/mean) of each respective parameter estimate over all noise realizations. Fifteen log-spaced SNR levels between 10 and 1000 are simulated. These scripts were generated using MATLAB (R2021b). 
