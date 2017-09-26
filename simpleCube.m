clear all
close all
clc
warning('off','all')

gravel(1).radius = 0.02;
gravel(1).coordinates = [0.25,0.5,0.5];
gravel(2).radius = 0.02;
gravel(2).coordinates = [0.5,0.25,0.5];
gravel(3).radius = 0.02;
gravel(3).coordinates = [0.5,0.5,0.25];

cube.corners = [0,0,0
                0,0,1
                0,1,0
                0,1,1
                1,0,0
                1,0,1
                1,1,0
                1,1,1];
cube.planes.xback.normal  = [-1, 0, 0];
cube.planes.xfront.normal = [ 1, 0, 0];
cube.planes.yback.normal  = [ 0,-1, 0];
cube.planes.yfront.normal = [ 0, 1, 0];
cube.planes.zback.normal  = [ 0, 0,-1];
cube.planes.zfront.normal = [ 0, 0, 1];
cube.shrinkRate = 0.000000000001;

% initiate speed of gravel
gravel(1).velocity = -1 + (1 --1)*rand(1,3);
gravel(1).velocity = gravel(1).velocity/norm(gravel(1).velocity);
% gravel(1).velocity = [1,0,0];
gravel(2).velocity = -1 + (1 --1)*rand(1,3);
gravel(2).velocity = gravel(2).velocity/norm(gravel(2).velocity);
% gravel(2).velocity = [0,-1,0];
gravel(3).velocity = -1 + (1 --1)*rand(1,3);
gravel(3).velocity = gravel(2).velocity/norm(gravel(2).velocity);
% gravel(3).velocity = [0,-1,0];
% gravel.speed = 0.1;

numberOfEvents = 20;
numberOfgravels = 3;

for i=1:numberOfEvents
    gravelCollisionTime = findMinimumgravelCollisionTime(gravel)
    for j=1:numberOfgravels
        gravel(j).storedCoordinates(i,:) = gravel(j).coordinates;
        % determine collition distance to closest wall and update position
        [collitionPosition,speed] = findDistance(cube,gravel(j));
        gravel(j).coordinates = gravel(j).coordinates + collitionPosition;
        gravel(j).velocity = speed.new;
    end
end

hold on
for i=1:numberOfgravels
    plot3(gravel(i).storedCoordinates(:,1),gravel(i).storedCoordinates(:,2),gravel(i).storedCoordinates(:,3))
    bubbleplot3(gravel(i).storedCoordinates(:,1),gravel(i).storedCoordinates(:,2),gravel(i).storedCoordinates(:,3),gravel(i).radius*ones(length(gravel(i).storedCoordinates),1));
end
axis([0 1 0 1 0 1])
axis square
grid on


function collisionTime = findMinimumgravelCollisionTime(gravel)
% create all gravel interactions
c = combnk(1:length(gravel),2);
collisionTime = inf*ones(length(c),1);
for i=1:length(c)
    % solving for time when gravels collide
    % http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment
    % https://en.wikipedia.org/wiki/Quadratic_equation
    wnull = gravel(c(i,1)).coordinates - gravel(c(i,2)).coordinates;
    
    p = 2*wnull*(gravel(c(i,1)).velocity - gravel(c(i,2)).velocity)'/norm(gravel(c(i,1)).velocity - gravel(c(i,2)).velocity)^2;
    q = (wnull*wnull' - (gravel(c(i,1)).radius + gravel(c(i,2)).radius)^2)/norm(gravel(c(i,1)).velocity - gravel(c(i,2)).velocity)^2;
    
    t.plus = -p/2 + sqrt((p/2)^2 - q);
    t.minus = -p/2 - sqrt((p/2)^2 - q);
    
    % only solutions that are real valued and positive
    if isreal(t.plus) && isnan(t.plus)~=1 && isinf(t.plus)~=1 && t.plus > 0
        if t.minus < 0
            collisionTime(i) = t.plus;
        else
            collisionTime(i) = t.minus;
        end
    end
end
collisionTime = min(collisionTime);
end


function [collitionPosition,speed] = findDistance(cube,gravel)
lnull = gravel.coordinates;
l = gravel.velocity;

k = 0;
for i=1:6
    % check which plan the gravel will intersect with
    % https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    if i==1
        pnull = cube.corners(1,:);
        n = cube.planes.xback.normal;
        cosangle = gravel.velocity*cube.planes.xback.normal'/(norm(gravel.velocity)*norm(cube.planes.xback.normal));
        
        % determine time until impact between gravel and moving plane
        deltat(i) = (cube.corners(1,1) - gravel.coordinates(1) - gravel.radius)/(gravel.velocity(1) - cube.shrinkRate);
    elseif i==2
        pnull = cube.corners(5,:);
        n = cube.planes.xfront.normal;
        cosangle = gravel.velocity*cube.planes.xfront.normal'/(norm(gravel.velocity)*norm(cube.planes.xfront.normal));
        
        % determine time until impact between gravel and moving plane
        deltat(i) = (cube.corners(5,1) - gravel.coordinates(1) - gravel.radius)/(gravel.velocity(1) - cube.shrinkRate);
    elseif i==3
        pnull = cube.corners(1,:);
        n = cube.planes.yback.normal;
        cosangle = gravel.velocity*cube.planes.yback.normal'/(norm(gravel.velocity)*norm(cube.planes.yback.normal));
        
        % determine time until impact between gravel and moving plane
        deltat(i) = (cube.corners(1,2) - gravel.coordinates(2) - gravel.radius)/(gravel.velocity(2) - cube.shrinkRate);
    elseif i==4
        pnull = cube.corners(3,:);
        n = cube.planes.yfront.normal;
        cosangle = gravel.velocity*cube.planes.yfront.normal'/(norm(gravel.velocity)*norm(cube.planes.yfront.normal));
        
        % determine time until impact between gravel and moving plane
        deltat(i) = (cube.corners(3,2) - gravel.coordinates(2) - gravel.radius)/(gravel.velocity(2) - cube.shrinkRate);
    elseif i==5
        pnull = cube.corners(1,:);
        n = cube.planes.zback.normal;
        cosangle = gravel.velocity*cube.planes.zback.normal'/(norm(gravel.velocity)*norm(cube.planes.zback.normal));
        
        % determine time until impact between gravel and moving plane
        deltat(i) = (cube.corners(1,3) - gravel.coordinates(3) - gravel.radius)/(gravel.velocity(3) - cube.shrinkRate);
    else
        pnull = cube.corners(2,:);
        n = cube.planes.zfront.normal;
        cosangle = gravel.velocity*cube.planes.zfront.normal'/(norm(gravel.velocity)*norm(cube.planes.zfront.normal));
        
        % determine time until impact between gravel and moving plane
        deltat(i) = (cube.corners(2,3) - gravel.coordinates(3) - gravel.radius)/(gravel.velocity(3) - cube.shrinkRate);
    end
    
    d = (pnull - lnull)*n'/(l*n');
    
    % adjust so that impact occures a gravel surface, not its center
    d = d - gravel.radius/cosangle;
    collitionVector = d*l + lnull - gravel.coordinates;
    if collitionVector*gravel.velocity' > 10*eps
        k = k + 1;
        collitionVectors(k,1:3) = collitionVector;
        
        % speed component perpendicular to the plane
        speed.perpendicular(k,:) = norm(gravel.velocity)*cosangle*n;
        speed.parallel(k,:) = gravel.velocity - speed.perpendicular(k,:);
        speed.new(k,:) = speed.parallel(k,:) + -speed.perpendicular(k,:);
    end
end

[collitionDistance,wall] = min(sqrt(sum(collitionVectors.^2,2)));
collitionPosition = collitionVectors(wall,:);
speed.new = speed.new(wall,:);
end