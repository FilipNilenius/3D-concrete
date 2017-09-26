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
gravel(1).velocity = [1,0,0];
gravel(2).velocity = -1 + (1 --1)*rand(1,3);
gravel(2).velocity = gravel(2).velocity/norm(gravel(2).velocity);
gravel(2).velocity = [0,-1,0];
gravel(3).velocity = -1 + (1 --1)*rand(1,3);
gravel(3).velocity = gravel(2).velocity/norm(gravel(2).velocity);
% gravel(3).velocity = [0,-1,0];
% gravel.speed = 0.1;

numberOfEvents = 2;
numberOfgravels = 2;

for i=1:numberOfEvents
    gravelCollisionTime = findMinimumgravelCollisionTime(gravel);
    [timeToWallCollision velocity] = findMinimumGravelWallCollisionTimee(cube,gravel);
    
    t = min([gravelCollisionTime timeToWallCollision])
    % update gravel coordinates at time of collision and speed(s) for
    % colliding gravels/gravel-wall
    for j=1:length(gravel)
        gravel(j).coordinates = gravel(j).coordinates + t*gravel(j).velocity;
        
        
        % store coordinates for plotting
        allCoordinates.gravel(j).event(i,:) = gravel(j).coordinates;
    end
end

hold on
for i=1:numberOfgravels
    plot3(allCoordinates.gravel(i).event(:,1),allCoordinates.gravel(i).event(:,2),allCoordinates.gravel(i).event(:,3))
    bubbleplot3(allCoordinates.gravel(i).event(:,1),allCoordinates.gravel(i).event(:,2),allCoordinates.gravel(i).event(:,3),gravel(i).radius*ones(length(allCoordinates.gravel(i).event),1));
end
axis([0 1 0 1 0 1])
axis square
% grid on
% hold off


function collisionTime = findMinimumgravelCollisionTime(gravel)
% create all gravel interaction combinations
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


function [timeToWallCollision collidingGravel velocity] = findMinimumGravelWallCollisionTimee(cube,gravel)
% all cominations between gravls and walls
% https://se.mathworks.com/matlabcentral/answers/98191-how-can-i-obtain-all-possible-combinations-of-given-vectors-in-matlab#answer_107541
[A,B] = meshgrid(1:length(gravel),1:6);
c=cat(2,A',B');
d=reshape(c,[],2);


xnull(1,:) = cube.corners(1,:);
xnull(2,:) = cube.corners(5,:);
xnull(3,:) = cube.corners(1,:);
xnull(4,:) = cube.corners(3,:);
xnull(5,:) = cube.corners(1,:);
xnull(6,:) = cube.corners(2,:);
n(1,:) = cube.planes.xback.normal;
n(2,:) = cube.planes.xfront.normal;
n(3,:) = cube.planes.yback.normal;
n(4,:) = cube.planes.yfront.normal;
n(5,:) = cube.planes.zback.normal;
n(6,:) = cube.planes.zfront.normal;

timeToWallCollision = inf*ones(length(d),1);
for i=1:length(d)
    timeToWallCollision(i) = (gravel(d(i,1)).radius + n(d(i,2),:)*(xnull(d(i,2),:) - gravel(d(i,1)).coordinates)')/(n(d(i,2),:)*gravel(d(i,1)).velocity');
end

% identify which gravel will collide with which wall
timeToWallCollision(timeToWallCollision<0) = inf;
[timeToWallCollision index] = min(timeToWallCollision);
collidingGravel = d(index,1);
wall = d(index,2);

% new speed after collision
% https://math.stackexchange.com/questions/1225494/component-of-a-vector-perpendicular-to-another-vector
velocity.old = gravel(collidingGravel).velocity;
velocity.perpendicular = gravel(collidingGravel).velocity*n(wall,:)'/(n(wall,:)*n(wall,:)')*n(wall,:);
velocity.parallel = velocity.old - velocity.perpendicular;
velocity.new = velocity.parallel + -velocity.perpendicular;
end