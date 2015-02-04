function [strain] = computePrincipalStrains(stressIn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% mean values in Voigt notation
strain.Voigt = strainIn;
strain.tensor = VoigtMap(strain.Voigt,'voigt2tensor');


% compute principal strains and sort them descendently
[strain.eigenvectors,strain.eigenvalues] = eig(strain.tensor);
[strain.eigenvalues,I] = sort(diag(strain.eigenvalues),'descend');
strain.eigenvectors = strain.eigenvectors(:,I);
end

