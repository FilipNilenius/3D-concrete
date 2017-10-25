function [t] = findTimeToWallCollision(xnull,x,v1,v2,r,n,d,t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(d)
    t(i) = (n(d(i,2),:)*(xnull(d(i,2),:) - x(d(i,1),:))' - r(d(i,1)))/(n(d(i,2),:)*(v1(d(i,1),:) - v2(d(i,2),:))');
end
end

