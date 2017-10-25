% QUADRATICEQUATION_SCRIPT   Generate MEX-function quadraticEquation_mex from
%  quadraticEquation.
% 
% Script generated from project 'quadraticEquation.prj' on 25-Oct-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;

%% Define argument types for entry-point 'quadraticEquation'.
ARGS = cell(1,1);
ARGS{1} = cell(6,1);
ARGS{1}{1} = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{6} = coder.typeof(0,[Inf  1],[1 0]);

%% Invoke MATLAB Coder.
codegen -config cfg quadraticEquation -args ARGS{1}
