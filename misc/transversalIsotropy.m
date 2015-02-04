function [stiffnessTensorVoigt] = transversalIsotropy(cement,ballast,ITZ,interfaceVoxel,meshProperties,iel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transversalIsotropy.m transforms the local transversal isotropic
% interface element to global coordinate system. The function returns
% element stiffness matrix in Voigt notation.
%
% Written by Filip Nilenius, 2013-08-29
% Last modified: 2013-09-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if interfaceVoxel.boundaryAggregate(iel) == 0
    E.ITZhat = 0.1/ITZ.YoungsModulus;
% else
%     E.ITZhat = 0;
% end
crack.one = interfaceVoxel.surfaceNormal(iel,:);


% base vecktors global system
e.one   = [1 0 0];
e.two   = [0 1 0];
e.three = [0 0 1];

% creates rotated and orthogonal coordinate system base on normal vector of
% elment
fony = crack.one + rand(1,3);
fony = fony/norm(fony);
crack.three = [crack.one(2)*fony(3)-crack.one(3)*fony(2) crack.one(3)*fony(1)-crack.one(1)*fony(3) crack.one(1)*fony(2)-crack.one(2)*fony(1)];
crack.three = crack.three/norm(crack.three);
crack.three = crack.three/norm(crack.three);
crack.two = [crack.three(2)*crack.one(3)-crack.three(3)*crack.one(2) crack.three(3)*crack.one(1)-crack.three(1)*crack.one(3) crack.three(1)*crack.one(2)-crack.three(2)*crack.one(1)];
crack.two = crack.two/norm(crack.two);


% creates transformation matrix
A11 = crack.one*e.one';
A12 = crack.one*e.two';
A13 = crack.one*e.three';
A21 = crack.two*e.one';
A22 = crack.two*e.two';
A23 = crack.two*e.three';
A31 = crack.three*e.one';
A32 = crack.three*e.two';
A33 = crack.three*e.three';

M = zeros(6,6); % for speed
M = [  A11^2     A12^2     A13^2       A12*A13         A11*A13         A11*A12
       A21^2     A22^2     A23^2       A22*A23         A21*A23         A21*A22
       A31^2     A32^2     A33^2       A32*A33         A31*A33         A31*A32
     2*A21*A31 2*A22*A32 2*A23*A33 A22*A33+A23*A32 A21*A33+A23*A31 A21*A32+A22*A31
     2*A11*A31 2*A12*A32 2*A13*A33 A12*A33+A13*A32 A13*A31+A33*A11 A11*A32+A12*A31
     2*A11*A21 2*A12*A22 2*A13*A23 A12*A23+A13*A22 A11*A23+A13*A21 A11*A22+A12*A21];


% Young's modulus. E.T = Reuss assumption, E.L = Voigt assumption

E.T = meshProperties.dx^3/(interfaceVoxel.volume.ballast(iel)/ballast.YoungsModulus + interfaceVoxel.volume.cement(iel)/cement.YoungsModulus + interfaceVoxel.area.ITZ(iel)*E.ITZhat);
E.L = interfaceVoxel.volume.ballast(iel)/meshProperties.dx^3*ballast.YoungsModulus + interfaceVoxel.volume.cement(iel)/meshProperties.dx^3*cement.YoungsModulus;% + t*interfaceVoxel.area.ITZ(iel)/meshProperties.dx^3*E.ITZ;


v.p  = cement.PoissonsRatio;
v.xp = cement.PoissonsRatio;
G.xp = E.T/(1+v.xp)/2;


% setup compliance matrix
S11 = 1/E.T;
S12 = -v.p/E.L;
S13 = -v.xp/E.L;
S21 = S12;
S22 = 1/E.L;
S23 = -v.xp/E.L;
S31 = S13;
S32 = S23;
S33 = S22;
S44 = 1/(E.L/(1+v.p)/2);
S55 = 1/G.xp;
S66 = 1/G.xp;

S = [S11 S12 S13   0   0   0
     S21 S22 S23   0   0   0
     S31 S32 S33   0   0   0
       0   0   0 S44   0   0
       0   0   0   0 S55   0
       0   0   0   0   0 S66];

C = inv(S);

% element stiffness matrix in global coordinate system
stiffnessTensorVoigt = M'*C*M;
end

