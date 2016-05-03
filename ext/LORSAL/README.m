% This set of files contains the Matlab code
% for the Bayesian supervised segmentation algorithm intoduced in the
% papers
%
%  J. Li, J. Bioucas-Dias, and Antonio Plaza, "Hyperspectral segmentation  
%  with active learning," submitted to IEEE TGRS,  2010 (link).
%
%  J. Li, J. Bioucas-Dias, and Antonio Plaza, "Semi-supervised Hyperspectral 
%  classification using active label selection," in SPIE Europe Remote 
%  Sensing, Vol.7477, 2009.
% 
%  J. Bioucas-Dias and M. Figueiredo,  "Logistic regression via variable 
%  splitting and augmented lagrangian tools"  Tech. Rep., Insituto Superior 
%  Tecnico, TULisbon, 2009
%
%
% 
% Files:
%                       readme.m   -  this file
%    demo_LORSAL_AL_MLL_AVIRIS.m   -  Matlab demo over real AVIRIS Indiana Pines dataset
%                LORSAL_MLL_AS.m   -  Matlab function file integrating
%                                     LORSAL and MLL spatial prior, and
%                                     adopting active sampling methods%                                   
%                       LORSAL.m   -  Matlab function file for LORSAL algorithm
%        train_test_random_new.m   -  Matlab function to randomly  select
%                                     training and testing set
%                      mlldata.m   -  Matlab functon to generate simulated
%                                     image following a MLL distribution.
%                    mlogistic.m   -  Matlab function to compute the
%                                      multinomial distributions (one per sample)
%            mlr_probabilities.m   -  Matlab function to compute the compute
%                                     class probalilities
%                    calcError.m   -  Matlab function to compute the
%                                     accuracies
%                     ToMatrix.m   -  Matlab function to reshape a vector
%                                     to a matrix 
%                     ToVector.m   -  Matlab function to reshape a matrix
%                                     to a vector
%  
%   AVIRIS_Indiana_16class.mat,imgreal.mat - Matlab data file with AVIRIS
%                                            Indiana Pines dataset
%
%
%
%             GCmex1.2.tar.gz   -  Matlab Wrapper for Graph Cuts
%                                  http://www.wisdom.weizmann.ac.il/~bagon/matlab.html
%
%%                 reference.bib   -  References. If you use this demo, please
%%                                  cite these references.
%
%  Any suggestions and comments are appreciated, and please send them to the
%  authors: jun@lx.it.pt. bioucas@lx.it.pt
%  August, 2010
%
%

