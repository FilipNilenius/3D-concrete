clear all
close all
clc
warning('off','all')

gravelSet(1) = gravel;
gravelSet(2) = gravel;
gravelSet(3) = gravel;

gravelSet(1).radius = 0.04;
gravelSet(1).coordinates = [0.26,0.5,0.5];
gravelSet(2).radius = 0.04;
gravelSet(2).coordinates = [0.5,0.25,0.5];
gravelSet(3).radius = 0.04;
gravelSet(3).coordinates = [0.25,0.5,0.25];

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

% initiate speed of gravelSet
gravelSet(1).velocity = -1 + (1 --1)*rand(1,3);
gravelSet(1).velocity = gravelSet(1).velocity/norm(gravelSet(1).velocity);
% gravelSet(1).velocity = [1,0,0];
gravelSet(2).velocity = -1 + (1 --1)*rand(1,3);
gravelSet(2).velocity = gravelSet(2).velocity/norm(gravelSet(2).velocity);
% gravelSet(2).velocity = [0,1,0];
gravelSet(3).velocity = -1 + (1 --1)*rand(1,3);
gravelSet(3).velocity = gravelSet(3).velocity/norm(gravelSet(2).velocity);
% gravelSet(3).velocity = [0,-1,0];
% gravelSet.speed = 0.1;

numberOfEvents = 60;
numberOfgravelSets = 3;


for j=1:length(gravelSet)
    allCoordinates.gravelSet(j).event(1,:) = gravelSet(j).coordinates;    
end

for i=1:numberOfEvents
    [gravelSetCollisionTime collidingGravels velocities] = findMinimumgravelSetCollisionTime(gravelSet);
    [timeToWallCollision collidingGravel velocity] = findMinimumgravelSetWallCollisionTimee(cube,gravelSet);
    
    
    gravelSetCollisionTime;
    timeToWallCollision;
    timeToNextEvent = min([gravelSetCollisionTime timeToWallCollision]);
    % update gravelSet coordinates at time of collision
    for j=1:length(gravelSet)
        gravelSet(j).coordinates = gravelSet(j).coordinates + timeToNextEvent*gravelSet(j).velocity;
        
        % store coordinates for plotting
        allCoordinates.gravelSet(j).event(i+1,:) = gravelSet(j).coordinates;
    end
    
    % update velocity after collision
    if timeToWallCollision < gravelSetCollisionTime % if next event is wall collision
        gravelSet(collidingGravel).velocity = velocity.new;
    else
        gravelSet(collidingGravels(1)).velocity = velocities.one.new;
        gravelSet(collidingGravels(2)).velocity = velocities.two.new;
    end
end

hold on
for i=1:numberOfgravelSets
    plot3(allCoordinates.gravelSet(i).event(:,1),allCoordinates.gravelSet(i).event(:,2),allCoordinates.gravelSet(i).event(:,3))
    bubbleplot3(allCoordinates.gravelSet(i).event(:,1),allCoordinates.gravelSet(i).event(:,2),allCoordinates.gravelSet(i).event(:,3),gravelSet(i).radius*ones(length(allCoordinates.gravelSet(i).event),1));
end
axis([0 1 0 1 0 1 0 1])
axis square
% grid on
% hold off


function [collisionTime collidingGravels velocity] = findMinimumgravelSetCollisionTime(gravelSet)
% create all gravelSet interaction combinations
c = combnk(1:length(gravelSet),2);
collisionTime = inf*ones(length(c),1);
for i=1:length(c)
    % solving for time when gravelSets collide
    % http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment
    % https://en.wikipedia.org/wiki/Quadratic_equation
    wnull = gravelSet(c(i,1)).coordinates - gravelSet(c(i,2)).coordinates;
    p = 2*wnull*(gravelSet(c(i,1)).velocity - gravelSet(c(i,2)).velocity)'/norm(gravelSet(c(i,1)).velocity - gravelSet(c(i,2)).velocity)^2;
    q = (wnull*wnull' - (gravelSet(c(i,1)).radius + gravelSet(c(i,2)).radius)^2)/norm(gravelSet(c(i,1)).velocity - gravelSet(c(i,2)).velocity)^2;
    
    t.plus = -p/2 + sqrt((p/2)^2 - q);
    t.minus = -p/2 - sqrt((p/2)^2 - q);
    
    % only solutions that are real valued and positive
    if isreal(t.plus) && isnan(t.plus)~=1 && isinf(t.plus)~=1 && t.plus > 100*eps
        if t.minus < 0
            collisionTime(i) = t.plus;
        else
            collisionTime(i) = t.minus;
        end
    end
end
[collisionTime index] = min(collisionTime);
collidingGravels = [c(index,1) c(index,2)];

% compute updated speeds after collision
Gnull = gravelSet(c(index,1)).velocity - gravelSet(c(index,2)).velocity;
n = gravelSet(c(index,1)).coordinates - gravelSet(c(index,2)).coordinates;
n = n/norm(n);
velocity.one.new = gravelSet(c(index,1)).velocity - (n*Gnull')*2*gravelSet(c(index,2)).mass/(gravelSet(c(index,1)).mass + gravelSet(c(index,2)).mass)*n;
velocity.two.new = gravelSet(c(index,2)).velocity + (n*Gnull')*2*gravelSet(c(index,1)).mass/(gravelSet(c(index,1)).mass + gravelSet(c(index,2)).mass)*n;
end


function [timeToWallCollision collidingGravel velocity] = findMinimumgravelSetWallCollisionTimee(cube,gravelSet)
% all cominations between gravls and walls
% https://se.mathworks.com/matlabcentral/answers/98191-how-can-i-obtain-all-possible-combinations-of-given-vectors-in-matlab#answer_107541
[A,B] = meshgrid(1:length(gravelSet),1:6);
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
    timeToWallCollision(i) = (-gravelSet(d(i,1)).radius + n(d(i,2),:)*(xnull(d(i,2),:) - gravelSet(d(i,1)).coordinates)')/(n(d(i,2),:)*gravelSet(d(i,1)).velocity');
end

% identify which gravelSet will collide with which wall
timeToWallCollision(timeToWallCollision<100*eps) = inf;
[timeToWallCollision index] = min(timeToWallCollision);
collidingGravel = d(index,1);
wall = d(index,2);

% new speed after collision
% https://math.stackexchange.com/questions/1225494/component-of-a-vector-perpendicular-to-another-vector
velocity.old = gravelSet(collidingGravel).velocity;
velocity.perpendicular = gravelSet(collidingGravel).velocity*n(wall,:)'/(n(wall,:)*n(wall,:)')*n(wall,:);
velocity.parallel = velocity.old - velocity.perpendicular;
velocity.new = velocity.parallel + -velocity.perpendicular;
end