function [t] = findTimeToWallCollision(xnull,x,v1,v2,r,n,d,t,combinationsRows)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for j=1:length(combinationsRows)
    i = combinationsRows(j);
    t(i) = 0.5*(n(d(i,2),:)*(xnull(d(i,2),:) - x(d(i,1),:))' - r(d(i,1)))/(n(d(i,2),:)*(v1(d(i,1),:) - v2(d(i,2),:))');
end

% identify which gravelSet will collide with which wall
t(t<0) = inf;
end

