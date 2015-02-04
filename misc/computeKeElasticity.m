function [Ke] = computeKeElasticity(ir,wp,detJ,B,D)

ngp = ir^3;
Ke = zeros(24,24);
BGauss = zeros(6,24);


for i=1:ngp
    BGauss = B(1+6*(i-1):6*i,:);
    Ke = Ke + BGauss'*D*BGauss*detJ*wp(i);
end