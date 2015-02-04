clear all
close all
% clc


A.t = 0.81;
A.c = 1.34;
B.t = 10450;
B.c = 2537;
epsilonNull = 1e-4;

gprim = 0.1;

kappa = linspace(epsilonNull,10*epsilonNull);
g.t = 1-(1-A.t)*epsilonNull./kappa - A.t*exp(-B.t*(kappa-epsilonNull));
g.c = 1-(1-A.c)*epsilonNull./kappa - A.c*exp(-B.c*(kappa-epsilonNull));

dgtdkappa = epsilonNull./kappa.^2 - A.t*epsilonNull./kappa.^2 + A.t*exp(B.t*epsilonNull)*B.t*exp(-B.t*kappa);
test = diff(g.t)./diff(kappa);

subplot(1,2,1)
hold on
plot(kappa,g.t)
plot(kappa,g.c)
hold off

subplot(1,2,2)
hold on
plot(kappa,dgtdkappa)
plot(kappa(1:end-1),test)