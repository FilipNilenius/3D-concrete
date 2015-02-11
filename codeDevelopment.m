%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I use this file for development only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

addpath([pwd,'/classFiles'])
addpath([pwd,'/misc'])
warning('off','all')



% % create object
SVE = SVEclass;

% % set properties
SVE.nx                = 40;
SVE.realizationNumber = 4;
SVE.Lbox              = 2;
SVE.aggFrac           = 0.42;

time.steps = 20;
time.stepsize = 0.1;

initialCondition = 1;


% % apply methods
SVE.setPath(); % uses current working directory
%SVE.setPath('C:\optional\path'); % if you want files to be saved elsewhere
% SVE.generateSVE(SVE.aggFrac,ballastRadii,gravelSieve,SVE.Lbox,domainFactor);
% SVE.meshSVE();
% SVE.writeTopology();
% SVE.LinElasticitySolver();
SVE.LinTransSolver(initialCondition,time);
SVE.TransPostProcessor();