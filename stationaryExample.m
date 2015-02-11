%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example which computes effective diffusvity
% tensor of an SVE. Comutations are based on
% stationary simmulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath([pwd,'/classFiles'])
addpath([pwd,'/misc'])
warning('off','all')



% % create object
SVE = SVEclass;

% % set properties
SVE.nx                = 40;
SVE.realizationNumber = 1;
SVE.Lbox              = 2;
SVE.aggFrac           = 0.3;

% % input for 'generateSVE()'
ballastRadii = [20 8 4 2]/2/10;   % Radius in [cm]. From http://www.sciencedirect.com/science/article/pii/S0168874X05001563
gravelSieve  = [.25 .25 .35 .15]; % Distribution in fraction. sum(gravelSieve) should be 1.0
aggFrac      = 0.30;              % Aggregate volume fraction
Lbox         = 6;                 % size of SVE [cm]
domainFactor = 1.5;               % ballast particles are distributed inside domainFactor*LBox


% % apply methods
SVE.setPath(); % uses current working directory
% SVE.setPath('C:\optional\path'); % if you want files to be saved elsewhere
SVE.generateSVE(SVE.aggFrac,ballastRadii,gravelSieve,SVE.Lbox,domainFactor);
SVE.meshSVE();
SVE.writeTopology();
SVE.computeEffectiveDiffusivtyTensor();
