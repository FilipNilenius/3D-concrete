function [out] = VoigtMap(in,flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Voigt2tensor
if strcmp('voigt2tensor',flag) == 1
    if size(in) == [6 6] % 4th order in in Voigt notation
        out = zeros(3,3,3,3);
        out(1,1,1,1) = in(1,1);
        out(2,2,1,1) = in(2,1);
        out(3,3,1,1) = in(3,1);
        out(2,3,1,1) = in(4,1);
        out(3,1,1,1) = in(5,1);
        out(1,2,1,1) = in(6,1);

        out(1,1,2,2) = in(1,2);
        out(2,2,2,2) = in(2,2);
        out(3,3,2,2) = in(3,2);
        out(2,3,2,2) = in(4,2);
        out(3,1,2,2) = in(5,2);
        out(1,2,2,2) = in(6,2);

        out(1,1,3,3) = in(1,3);
        out(2,2,3,3) = in(2,3);
        out(3,3,3,3) = in(3,3);
        out(2,3,3,3) = in(4,3);
        out(3,1,3,3) = in(5,3);
        out(1,2,3,3) = in(6,3);

        out(1,1,2,3) = in(1,4);
        out(2,2,2,3) = in(2,4);
        out(3,3,2,3) = in(3,4);
        out(2,3,2,3) = in(4,4);
        out(3,1,2,3) = in(5,4);
        out(1,2,2,3) = in(6,4);

        out(1,1,3,1) = in(1,5);
        out(2,2,3,1) = in(2,5);
        out(3,3,3,1) = in(3,5);
        out(2,3,3,1) = in(4,5);
        out(3,1,3,1) = in(5,5);
        out(1,2,3,1) = in(6,5);

        out(1,1,1,2) = in(1,6);
        out(2,2,1,2) = in(2,6);
        out(3,3,1,2) = in(3,6);
        out(2,3,1,2) = in(4,6);
        out(3,1,1,2) = in(5,6);
        out(1,2,1,2) = in(6,6);
    end
    if size(in) == [6 1] % 2th order in in Voigt notation
        out = zeros(3,3);
        out(1,1) = in(1,1);
        out(2,2) = in(2,1);
        out(3,3) = in(3,1);
        out(2,3) = in(4,1);
        out(3,2) = in(4,1);
        out(1,3) = in(5,1);
        out(3,1) = in(5,1);
        out(1,2) = in(6,1);
        out(2,1) = in(6,1);
    end
end

% tensor2Voigt
if strcmp('tensor2voigt',flag) == 1
    if size(in) == [3 3 3 3]
        out = zeros(6,6);
        out(1,1) = in(1,1,1,1);
        out(2,1) = in(2,2,1,1);
        out(3,1) = in(3,3,1,1);
        out(4,1) = in(2,3,1,1);
        out(5,1) = in(3,1,1,1);
        out(6,1) = in(1,2,1,1);

        out(1,2) = in(1,1,2,2);
        out(2,2) = in(2,2,2,2);
        out(3,2) = in(3,3,2,2);
        out(4,2) = in(2,3,2,2);
        out(5,2) = in(3,1,2,2);
        out(6,2) = in(1,2,2,2);

        out(1,3) = in(1,1,3,3);
        out(2,3) = in(2,2,3,3);
        out(3,3) = in(3,3,3,3);
        out(4,3) = in(2,3,3,3);
        out(5,3) = in(3,1,3,3);
        out(6,3) = in(1,2,3,3);

        out(1,4) = in(1,1,2,3);
        out(2,4) = in(2,2,2,3);
        out(3,4) = in(3,3,2,3);
        out(4,4) = in(2,3,2,3);
        out(5,4) = in(3,1,2,3);
        out(6,4) = in(1,2,2,3);

        out(1,5) = in(1,1,3,1);
        out(2,5) = in(2,2,3,1);
        out(3,5) = in(3,3,3,1);
        out(4,5) = in(2,3,3,1);
        out(5,5) = in(3,1,3,1);
        out(6,5) = in(1,2,3,1);

        out(1,6) = in(1,1,1,2);
        out(2,6) = in(2,2,1,2);
        out(3,6) = in(3,3,1,2);
        out(4,6) = in(2,3,1,2);
        out(5,6) = in(3,1,1,2);
        out(6,6) = in(1,2,1,2);
    end
end
end

