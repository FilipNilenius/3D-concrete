function [normalCrackStrain stiffnessTensorVoigt] = equivalentStrain(cement,strain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equivalentStrain.m computes
%
% Written by Filip Nilenius, 2013-08-29
% Last modified: 2013-09-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



P1 = strain.eigenvectors(:,1)*strain.eigenvectors(:,1)';
DedotdotP1 = zeros(3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
            DedotdotP1(i,j) = DedotdotP1(i,j) + cement.stiffnessTensor.full(i,j,k,l)*P1(k,l);
            end
        end
    end
end


thirdTerm.nominator = zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
            thirdTerm.nominator(i,j,k,l) = DedotdotP1(i,j)*DedotdotP1(k,l);
            end
        end
    end
end

thirdTerm.denominator = cement.lambda + 2*cement.shearModulus + cement.Ec;
thirdTerm.full = thirdTerm.nominator/thirdTerm.denominator;


Phat1 = zeros(3,3,3,3);
for kk=1:2
    test = zeros(3,3,3,3);
    p_sym = 1/2*(strain.eigenvectors(:,1)*strain.eigenvectors(:,kk+1)' + [strain.eigenvectors(:,1)*strain.eigenvectors(:,kk+1)']');
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    test(i,j,k,l) = p_sym(i,j)*p_sym(k,l);
                end
            end
        end
    end
    Phat1 = Phat1 + 2/(strain.eigenvalues(1)-strain.eigenvalues(kk+1))*test;
end

normalCrackStrain = 0;
for i=1:3
    for j=1:3
        normalCrackStrain = normalCrackStrain + DedotdotP1(i,j)*strain.tensor(i,j)/thirdTerm.denominator;
    end
end
Dec = cement.stiffnessTensor.full - 2*cement.shearModulus*normalCrackStrain*Phat1 - thirdTerm.full;
stiffnessTensorVoigt = VoigtMap(Dec,'tensor2voigt');

% stiffnessTensorVoigt
cement.stiffnessTensor.Voigt

test = VoigtMap(thirdTerm.full,'tensor2voigt')

test2 = VoigtMap(2*cement.shearModulus*normalCrackStrain*Phat1,'tensor2voigt')

end
