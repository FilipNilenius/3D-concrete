function [strain] = computeStresses(B,ed,ir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% strain and stress in Voit notation

strain.Voigt = zeros(6,ir^3);

for i=1:ir^3
    strain.Voigt(:,i) = B(1+6*(i-1):6*i,:)*ed;
end


strain.Voigt = [mean(strain.Voigt')]';
strain.tensor = VoigtMap(strain.Voigt,'voigt2tensor');


% compute principal strains and sort them descendently
[strain.eigenvectors,strain.eigenvalues] = eig(strain.tensor);
[strain.eigenvalues,I] = sort(diag(strain.eigenvalues),'descend');
strain.eigenvectors = strain.eigenvectors(:,I);
end

