clear all
clc

B = coder.typeof(0,[162,24],[1,0]);
wp = coder.typeof(0,[27,1],[1,0]);


% codegen computeStresses -args {zeros(6,6),B,zeros(24,1),zeros(1)}
codegen computeKeElasticity -args {zeros(1,1),wp,zeros(1,1),B,zeros(6,6)}
