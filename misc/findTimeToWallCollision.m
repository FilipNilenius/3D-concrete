function [timeToWallCollision] = findTimeToWallCollision(xnull,x,v1,v2,r,n,d)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

timeToWallCollision = inf*ones(length(d),1);
for i=1:length(d)
    timeToWallCollision(i) = (n(d(i,2),:)*(xnull(d(i,2),:) - x(d(i,1),:))' - r(d(i,1)))/(n(d(i,2),:)*(v1(d(i,1),:) - v2(d(i,2),:))');
end
end

