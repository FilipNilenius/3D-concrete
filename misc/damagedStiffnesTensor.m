function [stiffnessTensorVoigt,damage,equivalentStrain] = damagedStiffnesTensor(material,strain,damage,boundaryElement)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transversalIsotropy.m transforms the local transversal isotropic
% interface element to global coordinate system. The function returns
% element stiffness matrix in Voigt notation.
%
% Written by Filip Nilenius, 2013-08-29
% Last modified: 2013-09-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A.t = 0.81;
A.c = 1.34;
B.t = 10450;
B.c = 2537;



if norm(strain.eigenvalues) == 0 || boundaryElement == 1
    stiffnessTensorVoigt = material.stiffnessTensor.Voigt;
    equivalentStrain = 0;
    return
end

strainMcAuley = zeros(3,3);
for i = 1:3
    if strain.eigenvalues(i) > 0
        strainMcAuley = strainMcAuley + strain.eigenvalues(i)*strain.eigenvectors(:,i)*strain.eigenvectors(:,i)';
    end
end

equivalentStrain = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                equivalentStrain = equivalentStrain + strainMcAuley(i,j)*material.stiffnessTensor.full(i,j,k,l)*strainMcAuley(k,l);
            end
        end
    end
end
equivalentStrain = sqrt(equivalentStrain/material.YoungsModulus);


if equivalentStrain > material.crackingStrain && equivalentStrain > damage.kappa && damage.omega < 1.00
    damage.kappa = equivalentStrain;
    gprim = material.crackingStrain/damage.kappa^2 - A.t*material.crackingStrain/damage.kappa^2 + A.t*exp(B.t*material.crackingStrain)*B.t*exp(-B.t*damage.kappa);

    damage.omega = 1-(1-A.t)*material.crackingStrain/damage.kappa - A.t*exp(-B.t*(damage.kappa-material.crackingStrain));
    effectiveStress = VoigtMap(material.stiffnessTensor.Voigt*strain.Voigt,'voigt2tensor');

    effectiveStressCrossEffectiveStress = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    effectiveStressCrossEffectiveStress(i,j,k,l) = effectiveStress(i,j)*effectiveStress(k,l);
                end
            end
        end
    end
    stiffnessTensorVoigt = VoigtMap((1 - damage.omega)*material.stiffnessTensor.full - gprim/(material.YoungsModulus*equivalentStrain)*effectiveStressCrossEffectiveStress,'tensor2voigt');
else
    stiffnessTensorVoigt = (1 - damage.omega)*material.stiffnessTensor.Voigt;
end
end