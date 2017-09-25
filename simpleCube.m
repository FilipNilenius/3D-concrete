clear all
close all
clc

sphere(1).radius = 0.021;
sphere(1).coordinates = [0.25,0.5,0.5];
sphere(2).radius = 0.021;
sphere(2).coordinates = [0.5,0.25,0.5];
sphere(3).radius = 0.021;
sphere(3).coordinates = [0.5,0.5,0.25];

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

% initiate speed of sphere
sphere(1).velocity = -1 + (1 --1)*rand(1,3);
sphere(1).velocity = sphere(1).velocity/norm(sphere(1).velocity);
sphere(2).velocity = -1 + (1 --1)*rand(1,3);
sphere(2).velocity = sphere(2).velocity/norm(sphere(2).velocity);
sphere(3).velocity = -1 + (1 --1)*rand(1,3);
sphere(3).velocity = sphere(3).velocity/norm(sphere(3).velocity);
% sphere.speed = 0.1;

numberOfEvents = 3;
numberOfSpheres = 3;

for i=1:numberOfEvents
    for j=1:numberOfSpheres
        sphere(j).storedCoordiantes(i,:) = sphere(j).coordinates;
        % determine collition distance to closest wall and update position
        [collitionPosition,speed] = findDistance(cube,sphere(j));
        sphere(j).coordinates = sphere(j).coordinates + collitionPosition;
        sphere(j).velocity = speed.new;
    end
end

hold on
for i=1:numberOfSpheres
    plot3(sphere(i).storedCoordiantes(:,1),sphere(i).storedCoordiantes(:,2),sphere(i).storedCoordiantes(:,3))
    bubbleplot3(sphere(i).storedCoordiantes(:,1),sphere(i).storedCoordiantes(:,2),sphere(i).storedCoordiantes(:,3),sphere(i).radius*ones(length(sphere(i).storedCoordiantes),1));
end
axis([0 1 0 1 0 1])
axis square
grid on


function [collitionPosition,speed] = findDistance(cube,sphere)
lnull = sphere.coordinates;
l = sphere.velocity;

k = 0;
for i=1:6
    % check which plan the sphere will intersect with
    % https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    if i==1
        pnull = cube.corners(1,:);
        n = cube.planes.xback.normal;
        cosangle = sphere.velocity*cube.planes.xback.normal'/(norm(sphere.velocity)*norm(cube.planes.xback.normal));
        
        % determine time until impact between sphere and moving plane
        deltat(i) = (cube.corners(1,1) - sphere.coordinates(1) - sphere.radius)/(sphere.velocity(1) - cube.shrinkRate);
    elseif i==2
        pnull = cube.corners(5,:);
        n = cube.planes.xfront.normal;
        cosangle = sphere.velocity*cube.planes.xfront.normal'/(norm(sphere.velocity)*norm(cube.planes.xfront.normal));
        
        % determine time until impact between sphere and moving plane
        deltat(i) = (cube.corners(5,1) - sphere.coordinates(1) - sphere.radius)/(sphere.velocity(1) - cube.shrinkRate);
    elseif i==3
        pnull = cube.corners(1,:);
        n = cube.planes.yback.normal;
        cosangle = sphere.velocity*cube.planes.yback.normal'/(norm(sphere.velocity)*norm(cube.planes.yback.normal));
        
        % determine time until impact between sphere and moving plane
        deltat(i) = (cube.corners(1,2) - sphere.coordinates(2) - sphere.radius)/(sphere.velocity(2) - cube.shrinkRate);
    elseif i==4
        pnull = cube.corners(3,:);
        n = cube.planes.yfront.normal;
        cosangle = sphere.velocity*cube.planes.yfront.normal'/(norm(sphere.velocity)*norm(cube.planes.yfront.normal));
        
        % determine time until impact between sphere and moving plane
        deltat(i) = (cube.corners(3,2) - sphere.coordinates(2) - sphere.radius)/(sphere.velocity(2) - cube.shrinkRate);
    elseif i==5
        pnull = cube.corners(1,:);
        n = cube.planes.zback.normal;
        cosangle = sphere.velocity*cube.planes.zback.normal'/(norm(sphere.velocity)*norm(cube.planes.zback.normal));
        
        % determine time until impact between sphere and moving plane
        deltat(i) = (cube.corners(1,3) - sphere.coordinates(3) - sphere.radius)/(sphere.velocity(3) - cube.shrinkRate);
    else
        pnull = cube.corners(2,:);
        n = cube.planes.zfront.normal;
        cosangle = sphere.velocity*cube.planes.zfront.normal'/(norm(sphere.velocity)*norm(cube.planes.zfront.normal));
        
        % determine time until impact between sphere and moving plane
        deltat(i) = (cube.corners(2,3) - sphere.coordinates(3) - sphere.radius)/(sphere.velocity(3) - cube.shrinkRate);
    end
    
    d = (pnull - lnull)*n'/(l*n');
    
    % adjust so that impact occures a sphere surface, not its center
    d = d - sphere.radius/cosangle;
    collitionVector = d*l + lnull - sphere.coordinates;
    if collitionVector*sphere.velocity' > 10*eps
        k = k + 1;
        collitionVectors(k,1:3) = collitionVector;
        
        % speed component perpendicular to the plane
        speed.perpendicular(k,:) = norm(sphere.velocity)*cosangle*n;
        speed.parallel(k,:) = sphere.velocity - speed.perpendicular(k,:);
        speed.new(k,:) = speed.parallel(k,:) + -speed.perpendicular(k,:);
    end
end

[collitionDistance,wall] = min(sqrt(sum(collitionVectors.^2,2)));
collitionPosition = collitionVectors(wall,:);
speed.new = speed.new(wall,:);
end