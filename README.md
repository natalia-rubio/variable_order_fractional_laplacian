# variable_order_fractional_laplacian

Implementation of variable order fractional Laplacian technique described here (*link to VOFL document*) to compute gravitational potential.  This is joint work with Andrea Giusti, Marta D'Elia, Roberto Garappa, and Eric Darve.

To use, save all files to the same folder and open the *run_with_defaults* script.  Adjust the gravitational field (function handle s(|x|)), cell size, mass distribution radius, and domain size or leave defaults.  Toggle options to plot K (modified gravitational field), plot I (gravitational potential), plot orbital velocity, and compare to the explicit solution (constant s only) on and off as desired.

The error can be reduced by decreasing the cell size h and is O(h<sup>2</sup>).  Error can be alternatively decreased by increasing the number of quadrature points used to compute K and f (mass distribution/source) in *full_implementation*.


*note: MATLAB image processing and curve fitting toolboxes required*


**run_with_defaults** : script to set default parameters and run *full_implementation*

**full_implementation** : function to compute K, f, and I (gravitational potential).  Optional plotting of K, I, and orbital velocity.  Optional comparison to explicit solution (constant s only)

**compute_Kn0_3D** : function to compute K at cells not centered at the origin (one octant of domain)

**compute_K01_3D** : function to compute sphere-integral contribution to K at the origin-centered cell

**compute_K02_3D** : function to compute cube-integral contribution to K at the origin-centered cell

**ramp_func** : ramp function used for calculation of K at origin-centered cell

**compute_f** : function to compute f (one octant of domain)

**oct2full** : function to convert an octant to full domain

**convolve_3D_K** : function to perform convolution of K and f using discrete Fourier transform (DFT)

**K_ramp** : ramp function to suppress boundary values of K causing instability in DFT of K

**lgwt** : function to provide quadrature points and weights (written by Greg von Winckel - 02/25/2004)

**compute_FK_radial** : function to compute radial Fourier transform of K

**rad_fourier_int_exp** : function to compute Fourier transform of K for each frequency

**map_radial_to_3D** : function to map radial K and Fourier transform of K to 3D grid

**convolve_radial_FK** : function to perform convolution of K and f using discrete radial Fourier transform for K and discrete Fourier transform (DFT) for f

**explicit_sol** : compute analytical solution for constant s
