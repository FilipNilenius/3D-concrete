function [crackArea] = computeElementCrackArea(n,dx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes crack area from eigen strain
% http://en.wikipedia.org/wiki/Line-plane_intersection
%
% Written by Filip Nilenius, 2013-07-31
% Last edited on: 2013-09-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% line segment normal
normal = [1 0 0    % line 1
          0 1 0    % line 2
          1 0 0    % line 3
          0 1 0    % line 4
          1 0 0    % line 5
          0 1 0    % line 6
          1 0 0    % line 7
          0 1 0    % line 8
          0 0 1    % line 9
          0 0 1    % line 10
          0 0 1    % line 11
          0 0 1];  % line 12

% line segment starting point
lnoll = [0 0 0
         1 0 0
         0 1 0
         0 0 0
         0 0 1
         1 0 1
         0 1 1
         0 0 1
         0 0 0
         1 0 0
         1 1 0
         0 1 0]*dx;
     
% voxel's corner points
voxelCornerPoints = [0 0 0
                     1 0 0
                     1 1 0
                     0 1 0
                     0 0 1
                     1 0 1
                     1 1 1
                     0 1 1]*dx;

% center of voxel 
pnoll = [0.5 0.5 0.5]*dx; 


% determines intersection points between plane and voxel's line segments
intersectionPoints = zeros(10,3);
intersectionPointsCounter = 0;
for iline = 1:length(normal)
    d = ((pnoll - lnoll(iline,:))*n)/(normal(iline,:)*n);
    
    if d>=0 && d<=dx
        intersectionPointsCounter = intersectionPointsCounter + 1;
        intersectionPoints(intersectionPointsCounter,:) = d*normal(iline,:) + lnoll(iline,:);
    end
end
intersectionPoints(intersectionPointsCounter+1:end,:) = [];

% determine which side of the plane each voxel node is
sameSideCounter = 0;
oppositeSideCounter = 0;
for iPoint=1:length(voxelCornerPoints)
    if (voxelCornerPoints(iPoint,:) - pnoll)*n>0
        sameSideCounter = sameSideCounter + 1;
        point.sameSideofNormal(sameSideCounter,:) = voxelCornerPoints(iPoint,:);
    else
        oppositeSideCounter = oppositeSideCounter + 1;
        point.oppositeSideofNormal(oppositeSideCounter,:) = voxelCornerPoints(iPoint,:);
    end
end

% puts together each hull's point cloud
hull.belowNormal.pointCloud = [point.oppositeSideofNormal;intersectionPoints];
hull.aboveNormal.pointCloud = [point.sameSideofNormal;intersectionPoints];

% computes each hull's volume and area, respectively
[hull.belowNormal.volume,hull.belowNormal.area]=area3d(hull.belowNormal.pointCloud);
[hull.aboveNormal.volume,hull.belowNormal.area]=area3d(hull.aboveNormal.pointCloud);

crackArea = (hull.belowNormal.area + hull.belowNormal.area - 6*dx^2)/2;
end

