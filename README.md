Cross-Scale-Model-Integration-Examples
======================================

These scripts are provided to reproduce the examples given in the paper "Cross-scale 
integration of data and knowledge for predicting species ranges." All source code is
licensed under the GPL3, details are provided in the LICENSE file.

For convenience, a makefile is provided to demonstrate how to build the examples and
provide an automated way to generate all necessary files. The command *make all* will 
build both examples, while *make example1* and *make example2* will build each example 
individually. Note that these will likely take a long time and consume a lot of memory, 
so plan your usage accordingly.


Dependencies:
R version 3.x
JAGS

R packages:
sp
glm2
fields
rjags