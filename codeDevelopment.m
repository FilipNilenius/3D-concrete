%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Develop code using this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
clear all
close all

addpath([pwd,'/classFiles'])
addpath([pwd,'/misc'])
warning('off','all')



% % create object
SVE = SVEclass;

% % set properties
SVE.nx                = 70;
SVE.realizationNumber = 4;
SVE.Lbox              = 2;
SVE.aggFrac           = 0.42;
SVE.nGaussPoints      = 2;
SVE.strainIncrement   = 5e-8;
SVE.startLoadStep     = 2;
% SVE.endLoadStep       = 10;
SVE.H.grad(1)         = -1;
SVE.H.grad(2)         =  0;
SVE.H.grad(3)         =  0;



% % apply methods
SVE.setPath();
% SVE.meshSVE();
% SVE.getElasticityProperties();
% SVE.getBoundaryElements();
% SVE.writeTopology();
% SVE.LinElasticitySolver();
% SVE.printDamagedSVE(810);
% SVE.computeElementCrackArea(650:810);
SVE.homogenizedStress(1:810)
% SVE.LinStatSolver(759);
% SVE.diffTensorFunctionOfStrain(1:10:10)





