function [t] = quadraticEquation(x,v,r,c,t,gravelPairs)
% computes the coefficients p,q, for quadratic function
%
% x = coordinates for all gravel
% v = velocities for all gravel
% r = radii for all gravel
% c = combiantion matrix for all gravel

p = 0;
q = 0;
for j=1:length(gravelPairs)
    i = gravelPairs(j);
    % solving for time when gravelSets collide
    % http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment
    % https://en.wikipedia.org/wiki/Quadratic_equation
    wnull = x(c(i,1),:) - x(c(i,2),:);
    
    if norm(wnull) < 5*r(c(i,1)) % only check collision time within a 3 radii distance (for speed)
        p = 2*wnull*(v(c(i,1),:) - v(c(i,2),:))'/norm(v(c(i,1),:) - v(c(i,2),:))^2;
        q = (wnull*wnull' - (r(c(i,1)) + r(c(i,2)))^2)/norm(v(c(i,1),:) - v(c(i,2),:))^2;

        if (p/2)^2 - q > 0
            t(i) = min([-p/2 + sqrt((p/2)^2 - q) -p/2 - sqrt((p/2)^2 - q)]);
        end

        % only keep relevant times
        if t(i) < 0
            t(i) = inf;
        end
    end
end
end