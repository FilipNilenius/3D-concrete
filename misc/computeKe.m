function [Ke] = computeKe(ir,wp,detJ,B,diffusionTensor)

ngp = ir^3;
Ke = zeros(8,8);


for i=1:ngp
    Ke = Ke + B.gaussPoint{i}'*diffusionTensor*B.gaussPoint{i}*detJ.gaussPoint(i)*wp(i);
end