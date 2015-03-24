Cross-Scale-Model-Integration-Examples
======================================

These scripts are provided to reproduce the examples given in the paper "Cross-scale 
integration of data and knowledge for predicting species ranges." All source code is
licensed under the GPL3, details are provided in the LICENSE file.

For convenience, a makefile is provided for each example.

Example 1
---------
Building the model requires an installation of R, JAGS, and rjags.
For users familiar with Gnu Make, a makefile is provided as a guide.
To run the analyses, use ```make analysis```.
To produce the figures, use ```make figures```.
Additionally, global settings for controlling the analysis and the production of figures are contained in the file ```ex1_globals.r```.
There are a total of 6 steps to complete the analysis and build all figures.


Example 2
---------
Building the model requires an installation of R as well as a C++ compiler.
For users familiar with Gnu Make, a makefile is provided as a guide. 
To perform the compilation of the Metropolis sampler (and thus to check for errors in the build process that would prevent the model executing properly), use ```make model```.

