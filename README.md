Author: Steffen Limmer <steffen.limmer@tu-berlin.de>

Matlab code to learn structured nonlinear mmse estimators for signals distributed uniformly in generalized unit balls. The code was used to generate figure 4 and 5 presented in the paper specified below. Please cite the paper whenever you plan to use this code. To generate all results run the simulation scripts and 'plot_results.m' for Fig. 4, 'plot_pcoeffs.m' for Fig. 5.

`@inproceedings{LiSt16,
	author = {Limmer, S. and Stanczak, S.},
	title = {Towards optimal nonlinearities for sparse recovery using higher-order statistics},
	year = {2016},
	booktitle = {Proceedings of the IEEE International Workshop on Machine Learning For Signal Processing (MLSP)},
	month = {September 13--16},
	year = {2016},	
}`

[General Information]
- This software was tested on MATLAB Version R2014b using Manopt 2.0 (http://www.manopt.org) and CVX Version 1.22 (http://cvxr.com) under Linux 64bit
- the functions 'gamrnd.m', 'gen_vec.m', 'rndcheck.m' must be obtained from the Randomized Algorithms Control Toolbox (http://ract.sourceforge.net) (find them in '/ract/@ubase/private')
- the multinom* functions must be obtained from the GPLv3 package Specfun (http://octave.sourceforge.net/specfun/) 
- all credit for Manopt, CVX, RACT and Specfun goes to the original authors und significantly helped in the development of the code
- MATLAB (R) is a registered trademark of The MathWorks, Inc.
- Last updated: 25.07.2016

[Disclaimer]
The code is provided "as is", without warranty of any kind. In no event shall the authors or copyright holders be liable for any claim, damages or other liability. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.