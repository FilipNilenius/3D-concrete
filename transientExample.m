%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of transient simulation.
% Results are written to vtk-files.
% You need Paraview to visualize the results.
%
% The example takes approx. 5 min. to run.
%
% http://www.paraview.org/
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
SVE.nx                = 60;
SVE.realizationNumber = 1;
SVE.Lbox              = 2;
SVE.aggFrac           = 0.4;


% % input for 'generateSVE()'
ballastRadii = [20 8 4 2]/2/10;   % Radius in [cm]. From http://www.sciencedirect.com/science/article/pii/S0168874X05001563
gravelSieve  = [.25 .25 .35 .15]; % Distribution in fraction. sum(gravelSieve) should be 1.0
domainFactor = 1.5;               % ballast particles are distributed inside domainFactor*LBox


% % input for 'LinTransSolver()' method
time.steps = 10;
time.stepsize = 0.1;
initialCondition = 1;

% % set SVE boundary type
% % 1 = physical boundary where no aggregates cut the boundary surface
% % default value is 0 if not set by user
% SVE.boundary.x.back  = 1;
% SVE.boundary.x.front = 1;
% SVE.boundary.y.back  = 1;
% SVE.boundary.y.front = 1;
% SVE.boundary.z.back  = 1;
% SVE.boundary.z.front = 1;

% set convective coefficient at the three 'plus' faces of the SVE 
SVE.convCoeff.x.back = 100;
SVE.convCoeff.y.back = 100;
SVE.convCoeff.z.back = 100;


% % apply methods
SVE.setPath(); % uses current working directory
%SVE.setPath('C:\optional\path'); % if you want files to be saved elsewhere
SVE.generateSVE(SVE.aggFrac,ballastRadii,gravelSieve,SVE.Lbox,domainFactor);
SVE.meshSVE();
SVE.writeTopology();
SVE.LinTransSolver(initialCondition,time);
SVE.TransPostProcessor();
