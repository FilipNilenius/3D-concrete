classdef SVEclass < handle
    %----------------------------------------------------------------------
    % 3D model of heterogeneous concrete. Generates and anayses a
    % Statistical Volume Element (SVE) of mesoscale concrete. The models
    % includes methods to generate SVE realizations, meshing and
    % computations for both stationary and transient conditions.
    %
    % Model references: 
    % https://dx.doi.org/10.1007/s00466-014-0998-0
    % https://dx.doi.org/10.1007/s00466-014-1105-2
    %
    % Written by Filip Nilenius
    %----------------------------------------------------------------------
    properties
        % user defined
        realizationNumber
        aggFrac
        Lbox
        size
        nx
        transientName
        nGaussPoints = 2
        strainIncrement
        endLoadStep = 1000;
        startLoadStep
        residualTOL = 1e-3
        
        % set by constructor
        convCoeff
        ambHum
        H
        boundary
        
        % path to SVE realization on server
        path2Realization
    end
    methods(Access = public)
        function obj = SVEclass()
            % constructor
            obj.convCoeff.x.back  = inf;
            obj.convCoeff.x.front = inf;
            obj.convCoeff.y.back  = inf;
            obj.convCoeff.y.front = inf;
            obj.convCoeff.z.back  = inf;
            obj.convCoeff.z.front = inf;
            
            obj.ambHum.x.back  = 0;
            obj.ambHum.x.front = 0;
            obj.ambHum.y.back  = 0;
            obj.ambHum.y.front = 0;
            obj.ambHum.z.back  = 0;
            obj.ambHum.z.front = 0;
            
            obj.H.bar = 0;
            obj.H.grad(1) =  0;
            obj.H.grad(2) =  0;
            obj.H.grad(3) =  0;
            
            obj.boundary.x.back  = 0; % 1 if physical boundary (ie no gravel overlap)
            obj.boundary.x.front = 0;
            obj.boundary.y.back  = 0;
            obj.boundary.y.front = 0;
            obj.boundary.z.back  = 0;
            obj.boundary.z.front = 0;
            
            obj.transientName = 'defaultName';
        end
        function obj = setPath(obj,varargin)
            % setPath('string'):
            % Sets path to server. 'string' = 'bom' || 'cluster'.
            if nargin == 1
                Root = pwd;
            else
                Root = varargin{1};
            end
            
            obj.path2Realization = [Root,'/SVEs/Lbox_',num2str(obj.Lbox),'/aggFrac_',num2str(obj.aggFrac*100),'/Realization_',num2str(obj.realizationNumber),'/'];
            
            if ~exist(obj.path2Realization, 'dir')
                mkdir(obj.path2Realization);
            end
        end
        function generateSVE(obj,ballastRadii,gravelSieve,numberOfGravel)
            % generateSVE(obj,aggFrac,ballastRadii,gravelSieve,Lbox,domainFactor)
            %
            % generates 3D mesocale structure of concrete
            %    
            %   % example
            %   ballastRadii = [20 8 4 2]/2/10;   % Radius in [cm]. From http://www.sciencedirect.com/science/article/pii/S0168874X05001563
            %   gravelSieve  = [.25 .25 .35 .15]; % Distribution in fraction. sum(gravelSieve) should be 1.0
            %   aggFrac      = 0.30;              % Aggregate volume fraction
            %   Lbox         = 6;                 % size of SVE [cm]
            %   domainFactor = 1.5;               % ballast particles are distributed insidet domainFactor*LBox
            disp('---generating SVE---')
            
            % Break script if gravel sieve is wrongly defined
            if round(sum(gravelSieve),2) ~=1 || length(ballastRadii) ~= length(gravelSieve)
                error('gravel sieve wrongly defined')
            end
            
            % parameters
            cube.size = obj.Lbox;
            
            
            % create gravel set
            ballastMass = 4/3*pi*ballastRadii.^3;
            ratios = gravelSieve./ballastMass;
            numberOfGravelInSieve = round(ratios*numberOfGravel/sum(ratios));
            gravelSet(1:sum(numberOfGravelInSieve)) = gravel;
            cumSumGravelSieve = cumsum(numberOfGravelInSieve);
            
            
            for i=1:length(numberOfGravelInSieve)
                if i==1
                    range = 1:cumSumGravelSieve(i);
                else
                    range = cumSumGravelSieve(i-1)+1:cumSumGravelSieve(i);
                end
                [gravelSet(range).radius] = deal(ballastRadii(i));
            end
            gravelRadius = [gravelSet.radius];
            
            % define cube corners and normals
            cube.corners = [0,0,0
                            0,0,1
                            0,1,0
                            0,1,1
                            1,0,0
                            1,0,1
                            1,1,0
                            1,1,1]*cube.size;
            cube.planes.xback.normal  = [-1, 0, 0];
            cube.planes.xfront.normal = [ 1, 0, 0];
            cube.planes.yback.normal  = [ 0,-1, 0];
            cube.planes.yfront.normal = [ 0, 1, 0];
            cube.planes.zback.normal  = [ 0, 0,-1];
            cube.planes.zfront.normal = [ 0, 0, 1];
            
            
            % seed the random number generator based on the current time
            rng('shuffle')
            
            % distribute gravel inside cube
            k = 1;
            gravelCoordinates = zeros(length(gravelSet),3);
            while k <= length(gravelSet)
                bound.lower = 1.1*gravelSet(k).radius;
                bound.upper = cube.size - 1.1*gravelSet(k).radius;
                gravelCoordinates(k,:) = bound.lower + (bound.upper - bound.lower)*rand(1,3);
                if k >= 2
                    for i=1:k-1
                        gravelDistans = norm(gravelCoordinates(k,:) - gravelCoordinates(i,:));
                        sumOfRadii = gravelSet(k).radius + gravelSet(i).radius;
                        if gravelDistans < sumOfRadii
                            k = k - 1;
                            break
                        end
                    end
                end
                k = k + 1;
            end
            
            % Saves SVE data
            for i=1:length(gravelSet)
                gravelSet(i).coordinates = gravelCoordinates(i,:);
                gravelSet(i).radius = gravelRadius(i);
            end
            
            savefile = [obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'];
            save(savefile,'gravelSet','cube','gravelSieve');
        end
        function packSVE(obj,maxTime,varargin)
            % load SVE data
            load([obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'])
            
            % set shrink rate
            if isempty(varargin)==1
                cube.shrinkRate = -0.1;
            else
                cube.shrinkRate = 0;
            end
            
            % extract gravel coordinates
            gravelCoordinates = [reshape([gravelSet.coordinates],3,length(gravelSet))]';
            gravelRadius = [gravelSet.radius]';
            
            % create all gravelSet interaction combinations
            n = length(gravelSet);
            rows = (n - 1)*n/2;
            gravelCombinations = zeros(rows,2);
            
            counter = 0;
            for i=1:n-1
                for j=i+1:n
                    counter = counter + 1;
                    gravelCombinations(counter,1) = i;
                    gravelCombinations(counter,2) = j;
                end
            end
            clear n
            
            % create all gravel-wall interaction combinations
            [A,B] = meshgrid(1:length(gravelSet),1:6);
            c=cat(2,A',B');
            wallCombinations = reshape(c,[],2);
            
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
            v(1,:) = n(1,:)*0;
            v(2,:) = n(2,:)*cube.shrinkRate;
            v(3,:) = n(3,:)*0;
            v(4,:) = n(4,:)*cube.shrinkRate;
            v(5,:) = n(5,:)*0;
            v(6,:) = n(6,:)*cube.shrinkRate;
            
            % total gravel volume and assign random velocity
            totalVolume = 0;
            gravelVelocity = zeros(length(gravelSet),3);
            for j=1:length(gravelSet)
                totalVolume = totalVolume + gravelSet(j).mass;
                gravelSet(j).velocity = -1 + (1 --1)*rand(1,3);
                gravelSet(j).velocity = gravelSet(j).velocity/norm(gravelSet(j).velocity);
                gravelVelocity(j,:) = gravelSet(j).velocity;
            end
            
            % display volume fraction
            disp(['volume fraction before packing: ',num2str(totalVolume/(cube.size)^3)]);
            
            volumeFraction = 0;
            
            % pack aggregates
            disp('packing aggregates. please wait...')
            videoFile = VideoWriter('peaks.avi');
            open(videoFile);

            % define empty matrices
            gravelCollisionTimes = inf(length(gravelCombinations),1); % collision times for all gravel combinations
            wallCollisionTimes = inf(length(wallCombinations),1); % collision times for all gravel combinations'
            combinationRows.gravel = [1:length(gravelCombinations)]';
            combinationRows.wall = [1:length(wallCombinations)]';
            
            tic
            while volumeFraction < obj.aggFrac

                % finds min time of gravel-gravel collisions
                [gravelCollisionTimes] = quadraticEquation_mex(gravelCoordinates,gravelVelocity,...
                                                               gravelRadius,gravelCombinations,gravelCollisionTimes,combinationRows.gravel);
                [gravelCollisionTime gravelIndex] = min(gravelCollisionTimes);

                % finds time to next wall collision
                [wallCollisionTimes] = findTimeToWallCollision_mex(...
                                    xnull,...
                                    gravelCoordinates,...
                                    gravelVelocity,...
                                    v,...
                                    gravelRadius,...
                                    n,...
                                    wallCombinations,...
                                    wallCollisionTimes,...
                                    combinationRows.wall);
                [wallCollisionTime wallIndex] = min(wallCollisionTimes);

                % find time to next event
                timeToNextEvent = 0.8*min([gravelCollisionTime wallCollisionTime]); % 0.8 to avoid numerical errors

                % update gravelSet coordinates at time of collision
                gravelCoordinates = gravelCoordinates + timeToNextEvent*gravelVelocity;


                % new combinations matrices
                if wallCollisionTime < gravelCollisionTime
                    collidingGravel = wallCombinations(wallIndex,1);

                    one = find(gravelCombinations(:,1)==collidingGravel);
                    two = find(gravelCombinations(:,2)==collidingGravel);
                    combinationRows.gravel = [one' two']'; % indices in gravelCombination matrix where colliding gravel exists
                    combinationRows.wall = find(wallCombinations(:,1)==collidingGravel);
                else
                    collidingGravels = [gravelCombinations(gravelIndex,1) gravelCombinations(gravelIndex,2)];

                    subSet.gravel.one.left = find(gravelCombinations(:,1)==collidingGravels(1));
                    subSet.gravel.one.right = find(gravelCombinations(:,2)==collidingGravels(1));
                    subSet.gravel.one.full = [subSet.gravel.one.left' subSet.gravel.one.right']';
                    subSet.gravel.two.left = find(gravelCombinations(:,1)==collidingGravels(2));
                    subSet.gravel.two.right =  find(gravelCombinations(:,2)==collidingGravels(2)); % indices in gravelCombination matrix where colliding gravel exists
                    subSet.gravel.two.full = [subSet.gravel.two.left' subSet.gravel.two.right']';
                    subSet.wall.one = find(wallCombinations(:,1)==collidingGravels(1));
                    subSet.wall.two = find(wallCombinations(:,1)==collidingGravels(2));
                    combinationRows.wall = [subSet.wall.one' subSet.wall.two']';
                    combinationRows.gravel = [subSet.gravel.one.full' subSet.gravel.two.full']';    
                end

                % update time vectors
                gravelCollisionTimes = gravelCollisionTimes - timeToNextEvent;
                wallCollisionTimes = wallCollisionTimes - timeToNextEvent;
                gravelCollisionTimes(combinationRows.gravel) = inf;
                wallCollisionTimes(combinationRows.wall) = inf;


                % update cube corners (shrinking cube) at time of collision
                cube.corners(5,:) = cube.corners(5,:) + timeToNextEvent*cube.shrinkRate*cube.planes.xfront.normal;
                cube.corners(3,:) = cube.corners(3,:) + timeToNextEvent*cube.shrinkRate*cube.planes.yfront.normal;
                cube.corners(2,:) = cube.corners(2,:) + timeToNextEvent*cube.shrinkRate*cube.planes.zfront.normal;
                cube.size = cube.corners(5,1);
                xnull(1,:) = cube.corners(1,:);
                xnull(2,:) = cube.corners(5,:);
                xnull(3,:) = cube.corners(1,:);
                xnull(4,:) = cube.corners(3,:);
                xnull(5,:) = cube.corners(1,:);
                xnull(6,:) = cube.corners(2,:);


                % update gravel velocity after collision
                if wallCollisionTime < gravelCollisionTime % if next event is wall collision
                    % new speed after collision
                    % https://math.stackexchange.com/questions/1225494/component-of-a-vector-perpendicular-to-another-vector


                    collidingGravel = wallCombinations(wallIndex,1);
                    wall = wallCombinations(wallIndex,2);

                    velocity.old = gravelVelocity(collidingGravel,:);
                    velocity.perpendicular = gravelVelocity(collidingGravel,:)*n(wall,:)'/(n(wall,:)*n(wall,:)')*n(wall,:);
                    velocity.parallel = velocity.old - velocity.perpendicular;

                    % make sure that the perpendicular velocity is greater than the
                    % shrink rate of the cube. If not, make it 1.1 times shrink
                    % rate
                    if norm(velocity.perpendicular) < abs(cube.shrinkRate)
                        velocity.perpendicular = velocity.perpendicular/norm(velocity.perpendicular)*1.1*cube.shrinkRate;
                    end

                    gravelVelocity(collidingGravel,:) = velocity.parallel + -velocity.perpendicular;
                else % if next event is gravel-gravel collision
                    % compute updated speeds after collision
                    collidingGravels = [gravelCombinations(gravelIndex,1) gravelCombinations(gravelIndex,2)];
                    Gnull = gravelVelocity(gravelCombinations(gravelIndex,1),:) - gravelVelocity(gravelCombinations(gravelIndex,2),:);
                    normal = gravelCoordinates(gravelCombinations(gravelIndex,1),:) - gravelCoordinates(gravelCombinations(gravelIndex,2),:);
                    normal = normal/norm(normal);
                    gravelVelocity(collidingGravels(1),:) = gravelVelocity(gravelCombinations(gravelIndex,1),:) - (normal*Gnull')*2*gravelSet(gravelCombinations(gravelIndex,2)).mass/(gravelSet(gravelCombinations(gravelIndex,1)).mass + gravelSet(gravelCombinations(gravelIndex,2)).mass)*normal;
                    gravelVelocity(collidingGravels(2),:) = gravelVelocity(gravelCombinations(gravelIndex,2),:) + (normal*Gnull')*2*gravelSet(gravelCombinations(gravelIndex,1)).mass/(gravelSet(gravelCombinations(gravelIndex,1)).mass + gravelSet(gravelCombinations(gravelIndex,2)).mass)*normal;
                end

                % display volume fraction
                volumeFraction = totalVolume/(cube.size)^3;

%                 obj.makeMovie(videoFile,gravelSet,cube)
                if toc/(60*60) > maxTime
                    break
                end
            end
            % display volume fraction and SVe size
            disp(['volume fraction after packing: ',num2str(volumeFraction)]);
            disp(['SVE size after packing: ',num2str(cube.size)]);
            close(videoFile);
            close all

            % plot final state
            figure(2)
            bubbleplot3(gravelCoordinates(:,1),gravelCoordinates(:,2),gravelCoordinates(:,3),gravelRadius');
            xlim([0 cube.size])
            ylim([0 cube.size])
            zlim([0 cube.size])
            axis square
            
            % Saves SVE data
            for i=1:length(gravelSet)
                gravelSet(i).coordinates = gravelCoordinates(i,:);
                gravelSet(i).radius = gravelRadius(i);
            end
            
            savefile = [obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'];
            save(savefile,'gravelSet','cube','gravelSieve');
        end
        function meshSVE(obj)
            % meshSVE(nx):
            %   Mesh SVE with NEL = nx^3 and saves topology data to mat.file.
            %   nx = positive integer.
            disp('---meshing SVE---')
            
            tic
            load([obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'])
            obj.size = cube.size;
            dx = obj.size/(obj.nx-1);
            ny = obj.nx;
            nz = obj.nx;
            nel = obj.nx*ny*nz;
            
            % Voxel coordinates
            VoxelCoords = zeros(obj.nx*ny*nz,4);
            col3 = repmat(dx*[0:(obj.nx-1)],[obj.nx ny]);
            col4 = repmat(dx*[0:(obj.nx-1)],[obj.nx*ny 1]);
            VoxelCoords(:,2) = repmat(dx*[0:(obj.nx-1)],[1 ny*nz]);
            VoxelCoords(:,3) = col3(:)';
            VoxelCoords(:,4) = col4(:)';
            
            centroid = reshape([gravelSet.coordinates],3,length(gravelSet))';
            radius = [gravelSet.radius]';
            centro.x = centroid(:,1);
            centro.y = centroid(:,2);
            centro.z = centroid(:,3);

            centro.x = centro.x + radius;
            a = find(centro.x<0);
            centro.x = centro.x -2*radius;
            b = find(centro.x>obj.size);

            centro.y = centro.y + radius;
            c = find(centro.y<0);
            centro.y = centro.y -2*radius;
            d = find(centro.y>obj.size);

            centro.z = centro.z + radius;
            e = find(centro.z<0);
            centro.z = centro.z -2*radius;
            f = find(centro.z>obj.size);
            
            a = [a' b' c' d' e' f']';
            
            centroid(a,:) = [];
            radius(a) = [];
            
            % Determine which constituent (ballast-cement paste) each voxel is in
            ballastVoxels = zeros(length(VoxelCoords),1);
            for iballast=1:length(centroid(:,1))
                ccVector = bsxfun(@minus,VoxelCoords(:,2:end),centroid(iballast,:)); % c-c vector between current ballast and all other voxel elements
                normccVector = sqrt(sum(abs(ccVector).^2,2));
                diffNormRadii = normccVector - radius(iballast); % difference between norm and radius
                ballastVoxels = ballastVoxels + (diffNormRadii<0);
            end
            VoxelCoords(:,1) = ballastVoxels;


            % Creates topology matrix
            Edof = zeros(nel,9,'uint32');
            maxcoord = max(max(VoxelCoords(1:nel,2:4)));
            k = 0;
            for i=1:nel
                % diffusion
                Edof(i,1) = VoxelCoords(i,1);
                Edof(i,2) = i;
                Edof(i,3) = i  + 1;
                Edof(i,4) = obj.nx + i + 1;
                Edof(i,5) = obj.nx + i;
                Edof(i,6) = obj.nx*ny + i;
                Edof(i,7) = obj.nx*ny + i  + 1;
                Edof(i,8) = obj.nx*ny + obj.nx + i + 1;
                Edof(i,9) = obj.nx*ny + obj.nx + i;
                
                % Finds left-right-top boundary voxels to remove
                if VoxelCoords(i,2)>=maxcoord-dx/2 || VoxelCoords(i,3)>=maxcoord-dx/2 || VoxelCoords(i,4)>=maxcoord-dx/2
                    k = k + 1;
                    RemoveVoxel(1,k) = i;
                end
            end
            
            % Voxel center shifted to top rigt corner node
            NodeCoords = VoxelCoords;
            
            % Removes boundary voxels to ease global dof numbering
            Edof(RemoveVoxel,:) = [];

            meshProperties.ndof = max(max(Edof(:,2:end)));
            meshProperties.nel  = length(Edof(:,1));
            meshProperties.dx   = obj.size/(meshProperties.nel^(1/3));
            meshProperties.nx   = obj.nx;
            
            % Finds ITZ voxels
            [interfaceVoxel,Edof] = obj.findITZ(meshProperties,centroid,radius,NodeCoords,Edof);
            elementMaterial = Edof(:,1);
            
            % Determines boundary dofs
            plusFaceNodes  = find(max(NodeCoords(:,2:end)')>0.999*obj.size);
            minusFaceNodes = find(min(NodeCoords(:,2:end)')<0.001*obj.size);
            BoundaryNodes  = unique([plusFaceNodes minusFaceNodes]);
            
            % Determines 2D boundary elements along all faces
            [BEdof] = obj.boundaryEdof(Edof,meshProperties,NodeCoords,obj.size);
            
            % Generates sparse matrix structure % K = sparse(sparseVec.i,sparseVec.j,X,ndof,ndof)
            sparseVec.i = zeros(meshProperties.nel*8,8);
            sparseVec.j = zeros(meshProperties.nel*8,8);
            k=0;
            for i=1:meshProperties.nel
                % diffusion
                for j=1:8
                    k = k + 1;
                    sparseVec.i(k,:) = double(Edof(i,1+j))*ones(1,8);
                    sparseVec.j(k,:) = double(Edof(i,2:end));
                end
            end
            
            % diffusion
            sparseVec.i = reshape(sparseVec.i',meshProperties.nel*8*8,1);
            sparseVec.j = reshape(sparseVec.j',meshProperties.nel*8*8,1);
            
            % Write topology matrices to .mat-VoxelCoordss
            saveVoxelCoords = [obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'];
            save(saveVoxelCoords,'elementMaterial','Edof','BEdof','sparseVec','NodeCoords','BoundaryNodes','meshProperties','interfaceVoxel','-v7.3');
            
            % determine mesh properties pertaining to elasticity
            obj.getElasticityProperties();
            obj.getBoundaryElements();
            
            meshingTime = toc;
            if meshingTime > 60 && meshingTime < 60*60
                disp(['meshing done in ',num2str(meshingTime/60),' minutes'])
            elseif meshingTime >= 60*60
                disp(['meshing done in ',num2str(meshingTime/60/60),' hours'])
            else
                disp(['meshing done in ',num2str(meshingTime),' seconds'])
            end
        end
        function writeTopology(obj)
            % writeTopology():
            %   Write topology data to VTK-file.
            disp('---writes topology to file---')
            
            load ([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat']);
            fid = fopen([obj.path2Realization,'Topology_',num2str(meshProperties.nx),'_',num2str(obj.realizationNumber),'.vtk'], 'w');
            fprintf(fid,'# vtk DataFile Version 2.0\n');
            fprintf(fid,'Created on %s\n',datestr(now, 'yyyy-mm-dd'));
            fprintf(fid,'ASCII\n');
            fprintf(fid,'\n');
            fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
            fprintf(fid,'POINTS %d float\n',meshProperties.ndof.diffusion);
            fprintf(fid,'%f %f %f\n',NodeCoords(:,2:end)');
            fprintf(fid,'\n');
            fprintf(fid,'CELLS %d %d\n',meshProperties.nel,9*meshProperties.nel);
            fprintf(fid,'8 %d %d %d %d %d %d %d %d\n',Edof(:,2:end)'-1);
            fprintf(fid,'\n');
            fprintf(fid,'CELL_TYPES %d\n',meshProperties.nel);
            fprintf(fid,'%d\n',12*ones(meshProperties.nel,1));
            fprintf(fid,'\n');
            fprintf(fid,'CELL_DATA %d\n',meshProperties.nel);
            fprintf(fid,'SCALARS Material float 1\n');
            fprintf(fid,'LOOKUP_TABLE default\n');
            fprintf(fid,'%d\n',Edof(:,1));
            fclose(fid);
        end
        function [numAggFrac, numArea] = computeAggFracandITZarea(obj)
            % [numAggFrac numArea] = computeAggFracandITZarea():
            %   Returns aggregate fraction and total ITZ area.
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'interfaceVoxel','Edof','meshProperties')
            numAggFrac = 1 - (length(Edof(Edof(:,1)==0))*meshProperties.dx^3 + sum(interfaceVoxel.volume.cement))/obj.size^3;
            numArea = sum(interfaceVoxel.area.ITZ);
        end
        % stationary analysis
        function computeEffectiveDiffusivtyTensor(obj)
            % computeEffectiveDiffusivtyTensor():
            %   computes effective diffusivity tensor of SVE.
            disp('---solves stationary diffusion problem---')
            
            D = zeros(3,3);
            for i=1:3
                if i == 1
                    obj.H.grad(1) = -1;
                    obj.H.grad(2) =  0;
                    obj.H.grad(3) =  0;
                elseif i ==2
                    obj.H.grad(1) =  0;
                    obj.H.grad(2) = -1;
                    obj.H.grad(3) =  0;
                else
                    obj.H.grad(1) =  0;
                    obj.H.grad(2) =  0;
                    obj.H.grad(3) = -1;
                end
                D(:,i) = obj.LinStatSolver();
            end
            
            % try-catch if running without parallel toolbox
            try
                poolobj = gcp('nocreate');
                delete(poolobj);
            catch
            end
            dlmwrite([obj.path2Realization,'effective-diffusivity-tensor_nx=',num2str(obj.nx),'.txt'],D,'delimiter','\t')
        end
        function homoDiffRow = LinStatSolver(obj,varargin)
        % a = LinStatSolver(H,obj.sliceVector):
        %   Returns nodal vector 'a'.
        %   Solves the linear system Ka=f for a 3D mesostructure 
        %   of concrete. Dirichlet BC are implemented. Function is adapted for both
        %   full 3D SVE and 2D slices.
        
        
        % Load topology
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'interfaceVoxel')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'sparseVec');
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'BoundaryNodes')
        load([obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'],'cube')
        
        % needed for transparency for parfor
        elementMaterial = elementMaterial;
        interfaceVoxel = interfaceVoxel;
        obj.size = cube.size;
        
        crackDiffusivityTensor = cell(meshProperties.nel,1);
        crackDiffusivityTensor(:) = {zeros(3)};
        
        % load step solution if given as input
        nVarargs = length(varargin);
        if nVarargs == 1
            load([obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(varargin{1}),'.mat'])
            
            
            % compute crack widths in damaged elements
            
            [damegedElements test] = find(currentDamage > 0);
            elementCrackWidth = zeros(meshProperties.nel,1);
            
            
            for i = damegedElements'
                elementCrackWidth(i) = meshProperties.dx*maxPrincipalStrain.eigenvalue(i)*currentDamage(i);
            end
            
            max(elementCrackWidth)

            for iel = damegedElements'
                crackDiffusivityTensor{iel} = elementCrackArea(iel)*elementCrackWidth(iel)*cement.crackDiffusivity/(meshProperties.dx^3)*(eye(3) - maxPrincipalStrain.eigenvector(:,iel)*maxPrincipalStrain.eigenvector(:,iel)');
            end
        end
        
        
        % create material objects
        cement  = cementClass;
        ballast = ballastClass;
        ITZ     = ITZClass;
        
        % f = sparse(meshProperties.ndof,1);
        [ex_in,ey_in,ez_in] = obj.prepareKe(meshProperties.dx);
        [B,detJ,wp] = computeBmatrix(ex_in,ey_in,ez_in,obj.nGaussPoints);
        
        
        %X = zeros(meshProperties.nel*8,8);
        X = cell(meshProperties.nel,1);

        tic
        parfor iel=1:meshProperties.nel
            if elementMaterial(iel) == 2 % ITZ
                diffusionTensor = constitutiveModel(cement,ballast,ITZ,iel,interfaceVoxel,crackDiffusivityTensor{iel});
            elseif elementMaterial(iel) == 1 % Ballast
                diffusionTensor = ballast.diffusionCoefficient*eye(3);
            elseif elementMaterial(iel) == 0 % Cement
                diffusionTensor = cement.diffusionCoefficient*eye(3) + crackDiffusivityTensor{iel};
            end
            %X((iel*8-7):iel*8,:) = computeKe(2,wp,detJ,B,diffusionTensor);
            X{iel} = computeKe(2,wp,detJ,B,diffusionTensor);
        end
        
        assemblingK = toc;
        if assemblingK > 60 && assemblingK < 60*60
            disp(['assembling K in ',num2str(assemblingK/60),' minutes, finished on ', datestr(now)])
        elseif assemblingK >= 60*60
            disp(['assembling K in ',num2str(assemblingK/60/60),' hours, finished on ', datestr(now)])
        else
            disp(['assembling K in ',num2str(assemblingK),' seconds, finished on ', datestr(now)])
        end
        
        clear Edof;
        clear interfaceVoxel;
        
        X = cell2mat(X);
        X = reshape(X',meshProperties.nel*8*8,1);
        K.full = sparse(sparseVec.i,sparseVec.j,X,meshProperties.ndof.diffusion,meshProperties.ndof.diffusion);

        clear sparseVec;
        clear X;


        % Computes Dirichlet B.C. based on grad(v^M)
        a = zeros(meshProperties.ndof.diffusion,1);
        a(BoundaryNodes) = obj.H.grad(1)*(NodeCoords(BoundaryNodes,2) - 0.5*obj.size) +...
                           obj.H.grad(2)*(NodeCoords(BoundaryNodes,3) - 0.5*obj.size) +...
                           obj.H.grad(3)*(NodeCoords(BoundaryNodes,4) - 0.5*obj.size) +...
                           obj.H.bar;

        clear NodeCoords;

        % Condensates the system due to ballast dofs
        [dofs.cement.all, col] = find(diag(K.full)~=0);
        dofs.cement.prescribed = intersect(BoundaryNodes,dofs.cement.all);
        dofs.cement.free = setdiff(dofs.cement.all,dofs.cement.prescribed);
        K.fp = K.full(dofs.cement.free,dofs.cement.prescribed);
        K.ff = K.full(dofs.cement.free,dofs.cement.free);
        clear K.full;
        

        % Solves equation system using iterative solver
        tic
        a(dofs.cement.free) = minres(K.ff,-K.fp*a(dofs.cement.prescribed),[],70000,[],[],[]);
        solveK = toc;

        if solveK > 60 && solveK < 60*60
           disp(['solving Ka=f in ',num2str(solveK/60),' minutes'])
        elseif assemblingK >= 60*60
           disp(['solving Ka=f in ',num2str(solveK/60/60),' hours'])
        else
           disp(['solving Ka=f in ',num2str(solveK),' seconds'])
        end
            
        if nVarargs == 1
            homoDiffRow = obj.StatPostProcessor(a,currentDamage,maxPrincipalStrain,elementCrackArea);
        else
            homoDiffRow = obj.StatPostProcessor(a);
        end
        end
        function homoDiffRow = StatPostProcessor(obj,a,currentDamage,maxPrincipalStrain,elementCrackArea,varargin)
        % homoDiffRow = StatPostProcessor(a,obj.sliceVector,varargin):
        %   Computes homogenized diffusivity tensor components. If
        %   'varargin' exists, vector field is written to VTK.file (which
        %   requires a 'Topology_XX_X.vtk'-file.
        
        nVarargs = length(varargin);
        % Load topology
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'interfaceVoxel')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial')
        
        % needed for transparency for parfor
        Edof = Edof;
        elementMaterial = elementMaterial;
        interfaceVoxel = interfaceVoxel;
        
        % create material objects
        cement = cementClass;
        ballast = ballastClass;
        ITZ = ITZClass;
        
        crackDiffusivityTensor = cell(meshProperties.nel,1);
        crackDiffusivityTensor(:) = {zeros(3)};
        
        
        if nVarargs ~= 0
            % compute crack widths in damaged elements
            [damegedElements test] = find(currentDamage > 0);
            elementCrackWidth = zeros(meshProperties.nel,1);
        
            for i = damegedElements'
                elementCrackWidth(i) = meshProperties.dx*maxPrincipalStrain.eigenvalue(i)*currentDamage(i);
            end
            
            for iel = damegedElements'
                crackDiffusivityTensor{iel} = elementCrackArea(iel)*elementCrackWidth(iel)*cement.crackDiffusivity/(meshProperties.dx^3)*(eye(3) - maxPrincipalStrain.eigenvector(:,iel)*maxPrincipalStrain.eigenvector(:,iel)');
            end
        end
        
        % Computes flux in Gauss-points
        [ex_in,ey_in,ez_in] = obj.prepareKe(meshProperties.dx);
        [B] = computeBmatrix(ex_in,ey_in,ez_in,2);

        es_Store = zeros(meshProperties.nel,3);
        
        parfor iel=1:meshProperties.nel
            if elementMaterial(iel) == 2 % ITZ
                diffusionTensor = constitutiveModel(cement,ballast,ITZ,iel,interfaceVoxel,crackDiffusivityTensor{iel});
            elseif elementMaterial(iel) == 1 % Ballast
                diffusionTensor = ballast.diffusionCoefficient*eye(3);
            elseif elementMaterial(iel) == 0 % Cement
                diffusionTensor = cement.diffusionCoefficient*eye(3) + crackDiffusivityTensor{iel};
            end

            ed = a(Edof(iel,2:end))';

            elementFlux = zeros(3,1);
            for ii=1:(obj.nGaussPoints)^3
                elementFlux = elementFlux + -diffusionTensor*B.gaussPoint{ii}*ed';
            end

            es_Store(iel,:) = elementFlux/(obj.nGaussPoints)^3;
        end
        
        
        homoDiffRow = mean(es_Store)';
            
        
        % Write vector field to vtk-file
        if obj.H.bar(1) ~0
            writeVectorField(obj,es_Store,meshProperties.ndof,a,obj.path2Realization,obj.realizationNumber,obj.nx);
        end
        
        function writeVectorField(obj,es_Store,ndof,a,path2Realization,iRealization,nx)
        
        % Writes 3D structure data to .vtk-file
        fid = fopen([obj.path2Realization,'solution_I=0.vtk'], 'w');
        
        try
            s = fileread([obj.path2Realization,'Topology_',num2str(nx),'_',num2str(iRealization),'.vtk']);
        catch
            error('You first need to invoke the ''writeTopology()'' method.')
        end
        fprintf(fid,'%s',s);
        fprintf(fid, 'VECTORS vectors float\n');
        fprintf(fid, '%d %d %d\n',es_Store');
        fprintf(fid, '\n');
        fprintf(fid, 'POINT_DATA %d\n',meshProperties.ndof.diffusion);
        fprintf(fid, 'SCALARS stat float 1\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        fprintf(fid, '%d\n',a);
        fclose(fid);
        end
        end
        function diffTensorFunctionOfStrain(obj,strainVectorinput)
            
            % if alredy computed
            if exist([obj.path2Realization,'/diffTensorFunctionOfStrain.mat'])
                load([obj.path2Realization,'/diffTensorFunctionOfStrain.mat'])
            else
                diffTensor{1} = zeros(3,3);
                save([obj.path2Realization,'/diffTensorFunctionOfStrain.mat'],'diffTensor','-v7.3')
            end
            for strainIncrementVector = strainVectorinput
                load([obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(strainIncrementVector),'.mat'])
                
                for grad = 1:3
                    if grad == 1
                        obj.H.grad(1) = -1;
                        obj.H.grad(2) =  0;
                        obj.H.grad(3) =  0;
                    elseif grad == 2
                        obj.H.grad(1) =  0;
                        obj.H.grad(2) = -1;
                        obj.H.grad(3) =  0;
                    else
                        obj.H.grad(1) =  0;
                        obj.H.grad(2) =  0;
                        obj.H.grad(3) = -1;
                    end
                a = obj.LinStatSolver(strainIncrementVector);
                test(1:3,grad) = obj.StatPostProcessor(a,currentDamage,maxPrincipalStrain,elementCrackArea);
                end
                diffTensor{strainIncrementVector} = test;
                strainVector(strainIncrementVector) = totalStrain;
                save([obj.path2Realization,'/diffTensorFunctionOfStrain.mat'],'diffTensor','strainVector','-v7.3','-append')
            end
        end
        % transient analysis
        function LinTransSolver(obj,initialCondition,time)
            % LinTransSolver(initialCondition,time):
            %   Solves the linear system Cå+Ka=f for a 3D mesostrucutre 
            %   of concrete. Robin type BC are implemented.
            disp('---solves transient diffusion problem---')
            
            % Load topology
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'BEdof')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'interfaceVoxel')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'sparseVec');
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial')
            
            % needed for transparency for parfor
            elementMaterial = elementMaterial;
            interfaceVoxel = interfaceVoxel;
            
            % create material objects
            cement  = cementClass;
            ballast = ballastClass;
            ITZ     = ITZClass;
            

            Kce = meshProperties.dx^2*[ 1/9 3/54 1/36 3/54
                                       3/54  1/9 3/54 1/36
                                       1/36 3/54  1/9 3/54
                                       3/54 1/36 3/54  1/9];

            fce = meshProperties.dx^2*[1/4 1/4 1/4 1/4]';
            
            [ex_in,ey_in,ez_in] = obj.prepareKe(meshProperties.dx);
            [B,detJ,wp] = computeBmatrix(ex_in,ey_in,ez_in,2);
            Kee = obj.flw3i8e(ex_in,ey_in,ez_in,2,eye(3),2);
            Cee = obj.Ce_3D(ex_in,ey_in,ez_in,3,cement.storageCapacity);
            
            crackDiffusivityTensor = cell(meshProperties.nel,1);
            crackDiffusivityTensor(:) = {zeros(3)};

            % Assemble K and C
            K = cell(meshProperties.nel,1);
            C = cell(meshProperties.nel,1);
            
            tic
            parfor iel=1:meshProperties.nel
                if elementMaterial(iel) == 2 % ITZ
                    diffusionTensor = constitutiveModel(cement,ballast,ITZ,iel,interfaceVoxel,crackDiffusivityTensor{iel});
                    Ce = Cee;
                elseif elementMaterial(iel) == 1 % Ballast
                    diffusionTensor = ballast.diffusionCoefficient*eye(3);
                    Ce = 0*Cee;
                elseif elementMaterial(iel) == 0 % Cement
                    diffusionTensor = cement.diffusionCoefficient*eye(3) + crackDiffusivityTensor{iel};
                    Ce = Cee;
                end
                K{iel} = computeKe(2,wp,detJ,B,diffusionTensor);
                C{iel} = Ce;
            end
            poolobj = gcp('nocreate');
            delete(poolobj);
            
            K = cell2mat(K);
            K = reshape(K',meshProperties.nel*8*8,1);
            C = cell2mat(C);
            C = reshape(C',meshProperties.nel*8*8,1);
            
            K = sparse(sparseVec.i,sparseVec.j,K,meshProperties.ndof.diffusion,meshProperties.ndof.diffusion);
            C = sparse(sparseVec.i,sparseVec.j,C,meshProperties.ndof.diffusion,meshProperties.ndof.diffusion);
            f = sparse(meshProperties.ndof.diffusion,1);
            
            % Adds convective contribution to K and builds f
            for i=1:length(BEdof.x.back.Edof)
                if BEdof.x.back.Edof(i,5)~=2 && obj.convCoeff.x.back ~= inf % Excluding boundary dofs in ballast
                    K(BEdof.x.back.Edof(i,1:4),BEdof.x.back.Edof(i,1:4)) = K(BEdof.x.back.Edof(i,1:4),BEdof.x.back.Edof(i,1:4)) + obj.convCoeff.x.back*Kce;
                    f(BEdof.x.back.Edof(i,1:4)) = f(BEdof.x.back.Edof(i,1:4)) + obj.ambHum.x.back*obj.convCoeff.x.back*fce;
                end
                if BEdof.x.front.Edof(i,5)~=2 && obj.convCoeff.x.front ~= inf % Excluding boundary dofs in ballast
                    K(BEdof.x.front.Edof(i,1:4),BEdof.x.front.Edof(i,1:4)) = K(BEdof.x.front.Edof(i,1:4),BEdof.x.front.Edof(i,1:4)) + obj.convCoeff.x.front*Kce;
                    f(BEdof.x.front.Edof(i,1:4)) = f(BEdof.x.front.Edof(i,1:4)) + obj.ambHum.x.front*obj.convCoeff.x.front*fce;
                end
                if BEdof.y.back.Edof(i,5)~=2 && obj.convCoeff.y.back ~= inf % Excluding boundary dofs in ballast
                    K(BEdof.y.back.Edof(i,1:4),BEdof.y.back.Edof(i,1:4)) = K(BEdof.y.back.Edof(i,1:4),BEdof.y.back.Edof(i,1:4)) + obj.convCoeff.y.back*Kce;
                    f(BEdof.y.back.Edof(i,1:4)) = f(BEdof.y.back.Edof(i,1:4)) + obj.ambHum.y.back*obj.convCoeff.y.back*fce;
                end
                if BEdof.y.front.Edof(i,5)~=2 && obj.convCoeff.y.front ~= inf % Excluding boundary dofs in ballast
                    K(BEdof.y.front.Edof(i,1:4),BEdof.y.front.Edof(i,1:4)) = K(BEdof.y.front.Edof(i,1:4),BEdof.y.front.Edof(i,1:4)) + obj.convCoeff.y.front*Kce;
                    f(BEdof.y.front.Edof(i,1:4)) = f(BEdof.y.front.Edof(i,1:4)) + obj.ambHum.y.front*obj.convCoeff.y.front*fce;
                end
                if BEdof.z.back.Edof(i,5)~=2 && obj.convCoeff.z.back ~= inf % Excluding boundary dofs in ballast
                    K(BEdof.z.back.Edof(i,1:4),BEdof.z.back.Edof(i,1:4)) = K(BEdof.z.back.Edof(i,1:4),BEdof.z.back.Edof(i,1:4)) + obj.convCoeff.z.back*Kce;
                    f(BEdof.z.back.Edof(i,1:4)) = f(BEdof.z.back.Edof(i,1:4)) + obj.ambHum.z.back*obj.convCoeff.z.back*fce;
                end
                if BEdof.z.front.Edof(i,5)~=2 && obj.convCoeff.z.front ~= inf % Excluding boundary dofs in ballast
                    K(BEdof.z.front.Edof(i,1:4),BEdof.z.front.Edof(i,1:4)) = K(BEdof.z.front.Edof(i,1:4),BEdof.z.front.Edof(i,1:4)) + obj.convCoeff.z.front*Kce;
                    f(BEdof.z.front.Edof(i,1:4)) = f(BEdof.z.front.Edof(i,1:4)) + obj.ambHum.z.front*obj.convCoeff.z.front*fce;
                end
            end
            assem = toc;
            
            if assem > 60 && assem < 60*60
                disp(['assembling C and K in ',num2str(assem/60),' minutes'])
            elseif assem >= 60*60
                disp(['assembling C and K in ',num2str(assem/60/60),' hours'])
            else
                disp(['assembling C and K in ',num2str(assem),' seconds'])
            end
            
            % Condensates the system due to ballast dofs
            [CPdofs col] = find(diag(K)>0);
            K = K(CPdofs,CPdofs);
            C = C(CPdofs,CPdofs);
            f = f(CPdofs);
            a_int = initialCondition*ones(length(f),1);
            a_old = a_int;
            
            
            % Time iteration
            a_store = zeros(meshProperties.ndof.diffusion,time.steps);
            a_store(CPdofs,1) = initialCondition;
            for i=1:time.steps
                a_new = minres(C + time.stepsize*K,C*a_old + time.stepsize*f,[],70000,[],[],a_old);
                a_store(CPdofs,i+1) = a_new;
                a_old = a_new;
            end


            % Saves solution vector 'a' to file.
            savefile = [obj.path2Realization,obj.transientName,'.mat'];
            save(savefile,'a_store');

        end
        function TransPostProcessor(obj)
            % TransPostProcessor():
            %   Writes transient solution to VTK-file.
            
            disp('---post processing transient solution---')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'interfaceVoxel')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'BEdof')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial')
            load([obj.path2Realization,obj.transientName,'.mat']);
            
            % needed for transparency for parfor
            Edof = Edof;
            elementMaterial = elementMaterial;
            interfaceVoxel = interfaceVoxel;
            a_store = a_store;
            
            [ff,nts] = size(a_store); 
            
            
            % create material objects
            cement  = cementClass;
            ballast = ballastClass;
            ITZ     = ITZClass;
            
            crackDiffusivityTensor = cell(meshProperties.nel,1);
            crackDiffusivityTensor(:) = {zeros(3)};
            
            % Computes flux in Gauss-points
            [ex_in,ey_in,ez_in] = obj.prepareKe(meshProperties.dx);
            [B] = computeBmatrix(ex_in,ey_in,ez_in,2);


            for i=1:nts
                
                es_Store = zeros(meshProperties.nel,3);
                % Computes flux in Gauss-points
                parfor iel=1:meshProperties.nel
                    if elementMaterial(iel) == 2 % ITZ
                        diffusionTensor = constitutiveModel(cement,ballast,ITZ,iel,interfaceVoxel,crackDiffusivityTensor{iel});
                    elseif elementMaterial(iel) == 1 % Ballast
                        diffusionTensor = ballast.diffusionCoefficient*eye(3);
                    elseif elementMaterial(iel) == 0 % Cement
                        diffusionTensor = cement.diffusionCoefficient*eye(3) + crackDiffusivityTensor{iel};
                    end
                    
                    ed = a_store(Edof(iel,2:end),i)';
                    
                    elementFlux = zeros(3,1);
                    
                    for ii=1:(obj.nGaussPoints)^3
                        elementFlux = elementFlux + -diffusionTensor*B.gaussPoint{ii}*ed';
                    end

                    es_Store(iel,:) = elementFlux/(obj.nGaussPoints)^3;
                end
                
                % Writes 3D structure data to .vtk-file
                fid = fopen([obj.path2Realization,obj.transientName,'_',int2str(i-1),'.vtk'], 'w');
                s = fileread([obj.path2Realization,'Topology_',num2str(meshProperties.nx),'_',num2str(obj.realizationNumber),'.vtk']);
                fprintf(fid,'%s',s);
                fprintf(fid, 'VECTORS vectors float\n');
                fprintf(fid, '%d %d %d\n',es_Store');
                fprintf(fid, '\n');
                fprintf(fid, 'POINT_DATA %d\n',meshProperties.ndof.diffusion);
                fprintf(fid, 'SCALARS stat float 1\n');
                fprintf(fid, 'LOOKUP_TABLE default\n');
                fprintf(fid, '%d\n',a_store(:,i));
                fclose(fid);
                disp(['Number of time steps left: ',num2str(nts-i)])
            end
            poolobj = gcp('nocreate');
            delete(poolobj);
        end
        function [spatialDirection,slicedPlane] = plotTransientFront(obj)
            % [spatialDirection,slicedPlane] = plotTransientFront():
            %   Plots transient front for all time steps.
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
            load([obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'],'cube')
            load([obj.path2Realization,obj.transientName,'.mat']);
            
            
            [ff nts] = size(a_store);
            slicedPlane.concrete.mean = zeros(meshProperties.nx,1);
            slicedPlane.cement.std = zeros(meshProperties.nx,1);
            spatialDirection = linspace(0,cube.size,meshProperties.nx);
            
            for iTime=1:nts
                for i=0:meshProperties.nx-1
                    [nodeSlice foo] = find(abs(NodeCoords(:,2)-i*meshProperties.dx)<100*eps);
                    nodeSlice(nodeSlice==0) = [];
                    a.concrete = a_store(nodeSlice,iTime);
                    a.cement = a.concrete;
                    a.cement(a.cement==0) = [];
                    slicedPlane.concrete.mean(i+1) = mean(a.concrete);
                    slicedPlane.concrete.std(i+1)   = std(a.concrete);
                    slicedPlane.cement.mean(i+1) = mean(a.cement);
                    slicedPlane.cement.std(i+1)   = std(a.cement);
                end
                hold on
                figure (1)
                plot(spatialDirection,slicedPlane.concrete.mean)
                title('concrete (= cement + aggregates + ITZ)')
                xlabel('$x$-direction [cm]','Interpreter','latex')
                ylabel('Precentage of ambient concentration level')
                
                figure (2)
                title('cement and ITZ')
                plot(spatialDirection,slicedPlane.cement.mean)
                xlabel('$x$-direction [cm]','Interpreter','latex')
                ylabel('Precentage of ambient concentration level')
            end
        end
        % linear elasticity
        function LinElasticitySolver(obj)
        % LinElasticitySolver():
        %   solves linear elasticity problem with Dirichlet BC. returns
        %   solution vector 'a'.
        
        % Load topology
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'boundaryDofs');
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial');
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'interfaceVoxel');
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'EdofElasticity')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Nodedofs')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'BoundaryNodes')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'boundaryElement')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'BEdof')
        
        elementMaterial = elementMaterial;
        boundaryElement = boundaryElement;
        EdofElasticity = EdofElasticity;
        
        [ex,ey,ez] = obj.prepareKe(meshProperties.dx);
        [B,detJ,wp]=computeBmatrixElasticity(ex,ey,ez,obj.nGaussPoints);
        
        damage.omega = zeros(meshProperties.nel,1);
        damage.kappa = zeros(meshProperties.nel,1);
        
        cement = cementClass;
        ballast = ballastClass;
        ITZ = ITZClass;
        
        
        equivalentStrain = zeros(meshProperties.nel,1);
        maxPrincipalStrain.eigenvector = zeros(3,meshProperties.nel);
        maxPrincipalStrain.eigenvalue = zeros(1,meshProperties.nel);
        
        
        % compute stiffness tensors for all ITZ elements
        ITZelements = find(elementMaterial==2);
        stiffnessTensorITZ.Voigt = cell(meshProperties.nel,1);
        stiffnessTensorITZ.full = cell(meshProperties.nel,1);
        for i=1:length(ITZelements)
            stiffnessTensorITZ.Voigt{ITZelements(i)} = transversalIsotropy(cement,ballast,ITZ,interfaceVoxel,meshProperties,ITZelements(i));
            stiffnessTensorITZ.full{ITZelements(i)} = VoigtMap(transversalIsotropy(cement,ballast,ITZ,interfaceVoxel,meshProperties,ITZelements(i)),'voigt2tensor');
        end
        
        
        % apply boundary conditions
        a = zeros(meshProperties.ndof.elasticity,1);
        
        
        % hack
        apa1 = reshape(BEdof.z.back.Edof(:,1:4),1,4*length(BEdof.z.back.Edof));
        apa2 = reshape(BEdof.z.front.Edof(:,1:4),1,4*length(BEdof.z.back.Edof));
        BoundaryNodes = unique(sort([apa1 apa2]));
        % end hack
        
        dofs1 = Nodedofs(BoundaryNodes',1);
        cords.x = NodeCoords(BoundaryNodes',2);
        dofs2 = Nodedofs(BoundaryNodes',2);
        cords.y = NodeCoords(BoundaryNodes',3);
        dofs3 = Nodedofs(BoundaryNodes',3);
        cords.z = NodeCoords(BoundaryNodes',4);
        
        
        
        dofs.prescribed = [dofs1' dofs2' dofs3'];
        dofs.free = setdiff(1:meshProperties.ndof.elasticity,dofs.prescribed);
        
        
        
        
        % check if solutions from previous time steps exist
        if ~exist([obj.path2Realization,'/solution_',num2str(obj.nx)], 'dir') || obj.startLoadStep == 1
            obj.startLoadStep = 1;
            initiateLoadStep = 1;
            initiateRedTimeStep = 1;
            reducedStrainIncrement = obj.strainIncrement;
        else
            load([obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(obj.startLoadStep-1),'.mat'])
            damage.omega = currentDamage;
            damage.kappa = currentKappa;
            initiateLoadStep = obj.startLoadStep;
            initiateRedTimeStep = obj.startLoadStep;
        end
        
        
        
        loadStep = obj.startLoadStep;
        incrementSwitch = 1;
        while loadStep <= obj.endLoadStep
            disp(sprintf('\n'))
            disp(['load step ',num2str(loadStep),':'])
            
            strainedElements  = find(boundaryElement == 0);
            unDamagedElements = find(damage.omega == 0);
            
            
            undamagedElements = intersect(strainedElements,unDamagedElements);
            
            % update BC
            if loadStep == 1 || loadStep == 2
                totalStrain = obj.strainIncrement;
            elseif loadStep == 3
                % second load step at onset of cracking
                totalStrain = obj.strainIncrement + obj.strainIncrement*cement.crackingStrain/max(equivalentStrain);
            elseif incrementSwitch ==0
                totalStrain = totalStrain - obj.strainIncrement;
            else
                % if change in element damage is less than 1%%, scale def.
                % increment
                if max(damage.omega - oldDamage) < 0.01 && norm(damage.omega) > 0 && initiateLoadStep + 1 < loadStep
                    initiateLoadStep = loadStep;
                    reducedStrainIncrement = obj.strainIncrement;
                    totalStrain = totalStrain*cement.crackingStrain/max(equivalentStrain(undamagedElements));
                elseif max(damage.omega - oldDamage) > 0.60
                    reducedStrainIncrement = 0.1*reducedStrainIncrement;
                    if max(damage.omega - oldDamage) > 0.70
                        reducedStrainIncrement = 0.05*reducedStrainIncrement;
                    end
                    if max(damage.omega - oldDamage) > 0.80
                        reducedStrainIncrement = 0.01*reducedStrainIncrement;
                    end
                    if max(damage.omega - oldDamage) > 0.90
                        reducedStrainIncrement = 0.001*reducedStrainIncrement;
                    end
                    if reducedStrainIncrement < eps
                        reducedStrainIncrement = 0;
                    end
                    totalStrain = totalStrain + reducedStrainIncrement;
                else
                    totalStrain = totalStrain + reducedStrainIncrement;
                end
            end
            disp(['Total displacement: ',num2str(totalStrain*obj.size)])
            disp(['Total strain: ',num2str(totalStrain)])
            disp(['Displacement increment: ',num2str(reducedStrainIncrement)])
            
            % update BC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            u1 = -0.20*totalStrain*cords.x;
            u2 = -0.20*totalStrain*cords.y;
            u3 =       totalStrain*cords.z;
            
            a(dofs1) = u1;
            a(dofs2) = u2;
            a(dofs3) = u3;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            oldDamage = damage.omega;
            oldKappa = damage.kappa;
            residual = 1;
%             while norm(residual)>obj.residualTOL
            try
                for ii=1:3
                    % assemble K
                    K.full = sparse([],[],[],meshProperties.ndof.elasticity,meshProperties.ndof.elasticity);
                    
                    tic
                    for assemPart = 1:4
                        if assemPart == 1
                            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'firstSparseVecElasticity');
                            X = cell(length(firstSparseVecElasticity.first),1);
                            elements = firstSparseVecElasticity.first;
                            iVec = firstSparseVecElasticity.i;
                            jVec = firstSparseVecElasticity.j;
                            clear firstSparseVecElasticity
                        elseif assemPart ==2
                            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'secondSparseVecElasticity');
                            X = cell(length(secondSparseVecElasticity.second),1);
                            elements = secondSparseVecElasticity.second;
                            iVec = secondSparseVecElasticity.i;
                            jVec = secondSparseVecElasticity.j;
                            clear secondSparseVecElasticity
                        elseif assemPart ==3
                            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'thirdSparseVecElasticity');
                            X = cell(length(thirdSparseVecElasticity.third),1);
                            elements = thirdSparseVecElasticity.third;
                            iVec = thirdSparseVecElasticity.i;
                            jVec = thirdSparseVecElasticity.j;
                            clear thirdSparseVecElasticity
                        elseif assemPart == 4
                            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'fourthSparseVecElasticity');
                            X = cell(length(fourthSparseVecElasticity.fourth),1);
                            elements = fourthSparseVecElasticity.fourth;
                            iVec = fourthSparseVecElasticity.i;
                            jVec = fourthSparseVecElasticity.j;
                            clear fourthSparseVecElasticity
                        end

                        partEdof = EdofElasticity(elements,:);
                        partStiffnessTensorITZFull = stiffnessTensorITZ.full(elements);
                        partStiffnessTensorITZVoigt = stiffnessTensorITZ.Voigt(elements);
                        partElementMaterial = elementMaterial(elements);
                        partDamageOmega = damage.omega(elements);
                        partDamageKappa = damage.kappa(elements);
                        partBoundaryElement = boundaryElement(elements);
                        partEquivalentStrain = equivalentStrain(elements);

                        partMaxPrincipalEigenvector = cell(length(elements),1);
                        partMaxPrincipalStrain = cell(length(elements),1);

                        parfor iel=1:length(elements)
                            [elementStrain] = computeStresses(B.gaussPoint,a(partEdof(iel,:)),obj.nGaussPoints);
                            partMaxPrincipalEigenvector{iel} = elementStrain.eigenvectors(:,1);
                            partMaxPrincipalStrain{iel} = elementStrain.eigenvalues(1);

                            elementdamage.omega  = partDamageOmega(iel);
                            elementdamage.kappa  = partDamageKappa(iel);
                            if partElementMaterial(iel)==1 % Ballast
                                stiffnessTensorVoigt = ballast.stiffnessTensor.Voigt;

                            elseif partElementMaterial(iel)==2 % ITZ
                                material = ITZ;
                                material.stiffnessTensor.full = partStiffnessTensorITZFull{iel};
                                material.stiffnessTensor.Voigt = partStiffnessTensorITZVoigt{iel};
                                stiffnessTensorVoigt = partStiffnessTensorITZVoigt{iel};
                                [stiffnessTensorVoigt,elementdamage,partEquivalentStrain(iel)] = damagedStiffnesTensor(material,elementStrain,elementdamage,partBoundaryElement(iel));
                            elseif partElementMaterial(iel)==0 % Cement
                                stiffnessTensorVoigt = cement.stiffnessTensor.Voigt;
                                [stiffnessTensorVoigt,elementdamage,partEquivalentStrain(iel)] = damagedStiffnesTensor(cement,elementStrain,elementdamage,partBoundaryElement(iel));
                            end
                            partDamageOmega(iel) = elementdamage.omega;
                            partDamageKappa(iel) = elementdamage.kappa;
                            X{iel} = computeKeElasticity(obj.nGaussPoints,wp,detJ,B.gaussPoint,stiffnessTensorVoigt);
                        end

                        equivalentStrain(elements) = partEquivalentStrain;
                        damage.omega(elements) = partDamageOmega;
                        damage.kappa(elements) = partDamageKappa;

                        X = cell2mat(X);
                        X = reshape(X',length(elements)*8*3*8*3,1);
                        K.full = K.full + sparse(iVec,jVec,X,meshProperties.ndof.elasticity,meshProperties.ndof.elasticity);

                        partMaxPrincipalEigenvector = cell2mat(partMaxPrincipalEigenvector);
                        maxPrincipalStrain.eigenvector(:,elements) = reshape(partMaxPrincipalEigenvector,3,length(elements));
                        maxPrincipalStrain.eigenvalue(elements) = cell2mat(partMaxPrincipalStrain);

                    end
                    poolobj = gcp('nocreate');
                    delete(poolobj);
                    

                    clear X iVec jVec


                    assemblingK = toc;
                    if assemblingK > 60 && assemblingK < 60*60
                        disp(['assembling K in ',num2str(assemblingK/60),' minutes, finished on ', datestr(now)])
                    elseif assemblingK >= 60*60
                        disp(['assembling K in ',num2str(assemblingK/60/60),' hours, finished on ', datestr(now)])
                    else
                        disp(['assembling K in ',num2str(assemblingK),' seconds, finished on ', datestr(now)])
                    end

                    % create submatrices of K
                    K.fp = K.full(dofs.free,dofs.prescribed);
                    K.ff = K.full(dofs.free,dofs.free);
                    clear K.full;


                    % residual
                    residual = K.ff*a(dofs.free) + K.fp*a(dofs.prescribed);
                    disp(['||R||= ',num2str(norm(residual))])

                    % Solves for delta a_free(k), update a_free(k+1)
                    if norm(residual)>obj.residualTOL
                        tic
                        L1 = ichol(K.ff,struct('type','ict','droptol',1e-3));
                        deltaa = minres(K.ff,-residual,[],100000,L1,L1');

                        disp(['||deltaa||= ',num2str(norm(deltaa))])

                        solveK = toc;
                        a(dofs.free) = a(dofs.free) + deltaa;


                        if solveK > 60 && solveK < 60*60
                           disp(['solving Ka=f in ',num2str(solveK/60),' minutes'])
                        elseif solveK >= 60*60
                           disp(['solving Ka=f in ',num2str(solveK/60/60),' hours'])
                        else
                           disp(['solving Ka=f in ',num2str(solveK),' seconds'])
                        end

                        disp(['maximal damage: ',num2str(max(damage.omega))])
                        disp(['maximal delta damage: ',num2str(max(damage.omega - oldDamage))])
                        
                        if ii ~= 3
                            damage.omega = oldDamage;
                            damage.kappa = oldKappa;
                        end
                    end
                    clear L1 K;
                end
                
                % save current time step
                currentDamage = damage.omega;
                currentKappa = damage.kappa;

                % creat solution folder if needed
                if ~exist([obj.path2Realization,'/solution_',num2str(obj.nx)], 'dir')
                    mkdir([obj.path2Realization,'/solution_',num2str(obj.nx)]);
                end
                saveVoxelCoords = [obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(loadStep),'.mat'];
                save(saveVoxelCoords,'a','currentDamage','currentKappa','oldDamage','totalStrain','reducedStrainIncrement','maxPrincipalStrain','equivalentStrain','-v7.3');

                loadStep = loadStep + 1;
                incrementSwitch = 1;
            catch me
               incrementSwitch = 0;
            end
        end
        end
        function computeElementCrackArea(obj,chosenTimeSteps)
            load(['TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
            load(['TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial');
            
            for timeStep = chosenTimeSteps
                elementCrackArea = zeros(meshProperties.nel,1);
                load([obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(timeStep),'.mat'])
                for iel=1:meshProperties.nel
                    if currentDamage(iel) ~= 0 && elementMaterial(iel) == 0 % if cement or ITZ
                        crackNormal = maxPrincipalStrain.eigenvector(:,iel);
                        elementCrackArea(iel) = computeElementCrackArea(crackNormal,meshProperties.dx);
                    end
                end
                % save and append crack area to time step solution
                save([obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(timeStep),'.mat'],'elementCrackArea','-v7.3','-append');
            end
        end
        function printDamagedSVE(obj,loadSteps)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prints damaged SVE for timesteps in vector 'loadSteps' to
            % vtk-file for visualization in eg Paraview.
            %
            % Written by Filip Nilenius, 2013-08-15
            % Last edited on 2013-11-13
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial');
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
            
            % creat paraview folder if needed
            if ~exist([obj.path2Realization,'/paraview'], 'dir')
                mkdir([obj.path2Realization,'/paraview']);
            end
            
            for iLoadStep=loadSteps
                load([obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(iLoadStep)])
                
                % print solution to vtk-file
                fid = fopen([obj.path2Realization,'/paraview/Topology_LoadStep_',num2str(iLoadStep),'.vtk'], 'w');
%                 fid = fopen(['G:/Topology_LoadStep_',num2str(iLoadStep),'.vtk'], 'w');
                
                fprintf(fid,'# vtk DataFile Version 2.0\n');
                fprintf(fid,'Created on %s\n',datestr(now, 'yyyy-mm-dd'));
                fprintf(fid,'ASCII\n');
                fprintf(fid,'\n');
                fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
                fprintf(fid,'POINTS %d float\n',meshProperties.ndof.elasticity/3);
                fprintf(fid,'%f %f %f\n',NodeCoords(:,2:end)');
                fprintf(fid,'\n');
                fprintf(fid,'CELLS %d %d\n',meshProperties.nel,9*meshProperties.nel);
                fprintf(fid,'8 %d %d %d %d %d %d %d %d\n',Edof(:,2:end)'-1);
                fprintf(fid,'\n');
                fprintf(fid,'CELL_TYPES %d\n',meshProperties.nel);
                fprintf(fid,'%d\n',12*ones(meshProperties.nel,1));
                fprintf(fid,'\n');
                fprintf(fid,'CELL_DATA %d\n',meshProperties.nel);
                fprintf(fid,'SCALARS Material float 1\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                fprintf(fid,'%d\n',elementMaterial);
                fprintf(fid,'\n');
                fprintf(fid,'SCALARS damage float 1\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                fprintf(fid,'%d \n',currentDamage');
                fprintf(fid,'\n');
%                 fprintf(fid,'SCALARS strain_xx float 1\n');
%                 fprintf(fid,'LOOKUP_TABLE default\n');
%                 fprintf(fid,'%d\n',Strains(:,1)');
%                 fprintf(fid,'\n');
%                 fprintf(fid,'SCALARS strain_yy float 1\n');
%                 fprintf(fid,'LOOKUP_TABLE default\n');
%                 fprintf(fid,'%d\n',Strains(:,2)');
%                 fprintf(fid,'\n');
%                 fprintf(fid,'SCALARS strain_zz float 1\n');
%                 fprintf(fid,'LOOKUP_TABLE default\n');
%                 fprintf(fid,'%d\n',Strains(:,3)');
%                 fprintf(fid,'\n');
%                 fprintf(fid,'SCALARS stress_xx float 1\n');
%                 fprintf(fid,'LOOKUP_TABLE default\n');
%                 fprintf(fid,'%d\n',Stresses(:,1)');
%                 fprintf(fid,'\n');
%                 fprintf(fid,'SCALARS stress_yy float 1\n');
%                 fprintf(fid,'LOOKUP_TABLE default\n');
%                 fprintf(fid,'%d\n',Stresses(:,2)');
%                 fprintf(fid,'\n');
%                 fprintf(fid,'SCALARS stress_zz float 1\n');
%                 fprintf(fid,'LOOKUP_TABLE default\n');
%                 fprintf(fid,'%d\n',Stresses(:,3)');
%                 fprintf(fid,'\n');
                fprintf(fid,'SCALARS max_eigenstrain float 1\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                fprintf(fid,'%d\n', maxPrincipalStrain.eigenvalue');
                fprintf(fid,'\n');
                fprintf(fid,'SCALARS equivalent_strain float 1\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                fprintf(fid,'%d\n', equivalentStrain);
                fprintf(fid,'\n');
                fprintf(fid,'POINT_DATA %d\n',meshProperties.ndof.elasticity/3);
                fprintf(fid,'VECTORS deformation float\n');
                fprintf(fid,'%d %d %d\n',reshape(a,3,meshProperties.ndof.elasticity/3));
                fprintf(fid,'\n');
                fclose(fid);
            end
        end
        function homogenizedStress(obj,loadSteps)
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'elementMaterial');
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'interfaceVoxel')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'equivalentStrain')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'EdofElasticity')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'currentDamage')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'currentKappa')
            load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'boundaryElement')
            
            [ex,ey,ez] = obj.prepareKe(meshProperties.dx);
            [B,detJ,wp]=computeBmatrixElasticity(ex,ey,ez,obj.nGaussPoints);
            
            matlabpool open
            for iLoadStep=loadSteps
                load([obj.path2Realization,'/solution_',num2str(obj.nx),'/timeStepSolution_',num2str(obj.nx),'_',num2str(iLoadStep)])
                a = a;
                currentDamage = currentDamage;
                currentKappa = currentKappa;
                EdofElasticity = EdofElasticity;
                elementMaterial = elementMaterial;
                interfaceVoxel = interfaceVoxel;
                boundaryElement = boundaryElement;
                meshProperties = meshProperties;
                
                disp(['load step #: ',num2str(iLoadStep)])
                
                % if alredy computed
                if exist([obj.path2Realization,'/stressStrainRelation.mat'])
                    load([obj.path2Realization,'/stressStrainRelation.mat'])
                else
                    stressStrainRelation.stressTensorVoigt(1:6,1) = zeros(6,1);
                    stressStrainRelation.totalStrain(1) = 0;
                    save([obj.path2Realization,'/stressStrainRelation.mat'],'stressStrainRelation','-v7.3')
                end
                
                cement = cementClass;
                ballast = ballastClass;
                ITZ = ITZClass;
                
                
                
                stressTensorVoigt = zeros(6,meshProperties.nel);
                
                % compute stiffness tensors for all ITZ elements
                ITZelements = find(elementMaterial==2);
                stiffnessTensorITZ.Voigt = cell(meshProperties.nel,1);
                stiffnessTensorITZ.full = cell(meshProperties.nel,1);
                for i=1:length(ITZelements)
                    stiffnessTensorITZ.Voigt{ITZelements(i)} = transversalIsotropy(cement,ballast,ITZ,interfaceVoxel,meshProperties,ITZelements(i));
                    stiffnessTensorITZ.full{ITZelements(i)} = VoigtMap(transversalIsotropy(cement,ballast,ITZ,interfaceVoxel,meshProperties,ITZelements(i)),'voigt2tensor');
                end
                
                parfor iel=1:meshProperties.nel
                    if boundaryElement(iel) ~= 1 % not boundary elements
                        [elementStrain] = computeStresses(B.gaussPoint,a(EdofElasticity(iel,:)),obj.nGaussPoints);

                        elementdamage.omega  = currentDamage(iel);
                        elementdamage.kappa  = currentKappa(iel);
                        
                        if elementMaterial(iel)==1 % Ballast
                            stiffnessTensorVoigt = ballast.stiffnessTensor.Voigt;
                        elseif elementMaterial(iel)==2 % ITZ
                            material = ITZ;
                            material.stiffnessTensor.Voigt = stiffnessTensorITZ.Voigt{iel};
                            stiffnessTensorVoigt = stiffnessTensorITZ.Voigt{iel};
                            [stiffnessTensorVoigt,elementdamage,equivalentStrain(iel)] = damagedStiffnesTensor(material,elementStrain,elementdamage,boundaryElement(iel));
                        elseif elementMaterial(iel)==0 % Cement
                            stiffnessTensorVoigt = cement.stiffnessTensor.Voigt;
                            [stiffnessTensorVoigt,elementdamage,equivalentStrain(iel)] = damagedStiffnesTensor(cement,elementStrain,elementdamage,boundaryElement(iel));
                        end
                    stressTensorVoigt(:,iel) = stiffnessTensorVoigt*elementStrain.Voigt;    
                    end
                end
                
                stressStrainRelation.stressTensorVoigt(:,iLoadStep) = [mean(stressTensorVoigt')]';
                stressStrainRelation.totalStrain(iLoadStep) = totalStrain;
                
                save([obj.path2Realization,'/stressStrainRelation.mat'],'stressStrainRelation','-v7.3','-append');
            end
            matlabpool close
        end
        % update existing model with new functionallity
        function getElasticityProperties(obj)
        % getElasticityProperties():
        %   computes mesh data for elasticity given a mesh for
        %   diffusivity. Should only be used for existing mesh data that
        %   was done before elasticity was implemented.

        load([obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'])
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
        
        elementMaterial = Edof(:,1);
        nx = obj.nx;
        dx = obj.size/(nx-1);
        ny = nx;
        nz = nx;
        nel = nx*ny*nz;

        % Voxel coordinates
        VoxelCoords = zeros(nx*ny*nz,4);
        col3 = repmat(dx*[0:(nx-1)],[nx ny]);
        col4 = repmat(dx*[0:(nx-1)],[nx*ny 1]);
        VoxelCoords(:,2) = repmat(dx*[0:(nx-1)],[1 ny*nz]);
        VoxelCoords(:,3) = col3(:)';
        VoxelCoords(:,4) = col4(:)';


        % Creates topology matrix
        EdofElasticity = zeros(nel,8*3,'uint32');
        maxcoord = max(max(VoxelCoords(1:nel,2:4)));
        k = 0;
        for i=1:nel
            % elasticity
            % node 1
            EdofElasticity(i,1)  = 3*(i-1) + 1;
            EdofElasticity(i,2)  = 3*(i-1) + 2;
            EdofElasticity(i,3)  = 3*(i-1) + 3;

            % node 2
            EdofElasticity(i,4)  = 3*i + 1;
            EdofElasticity(i,5)  = 3*i + 2;
            EdofElasticity(i,6)  = 3*i + 3;

            % node 3
            EdofElasticity(i,7)  = nx*3 + 3*i + 1;
            EdofElasticity(i,8)  = nx*3 + 3*i + 2;
            EdofElasticity(i,9) = nx*3 + 3*i + 3;

            % node 4
            EdofElasticity(i,10) = nx*3 + 3*(i-1) + 1;
            EdofElasticity(i,11) = nx*3 + 3*(i-1) + 2;
            EdofElasticity(i,12) = nx*3 + 3*(i-1) + 3;

            % node 5
            EdofElasticity(i,13) = nx*ny*3 + 3*(i-1) + 1;
            EdofElasticity(i,14) = nx*ny*3 + 3*(i-1) + 2;
            EdofElasticity(i,15) = nx*ny*3 + 3*(i-1) + 3;

            % node 6
            EdofElasticity(i,16) = nx*ny*3 + 3*i + 1;
            EdofElasticity(i,17) = nx*ny*3 + 3*i + 2;
            EdofElasticity(i,18) = nx*ny*3 + 3*i + 3;

            % node 7
            EdofElasticity(i,19) = nx*ny*3 + nx*3 + 3*i + 1;
            EdofElasticity(i,20) = nx*ny*3 + nx*3 + 3*i + 2;
            EdofElasticity(i,21) = nx*ny*3 + nx*3 + 3*i + 3;

            % node 8
            EdofElasticity(i,22) = nx*ny*3 + nx*3 + 3*(i-1) + 1;
            EdofElasticity(i,23) = nx*ny*3 + nx*3 + 3*(i-1) + 2;
            EdofElasticity(i,24) = nx*ny*3 + nx*3 + 3*(i-1) + 3;

            % Finds left-right-top boundary voxels to remove
            if VoxelCoords(i,2)>=maxcoord-dx/2 || VoxelCoords(i,3)>=maxcoord-dx/2 || VoxelCoords(i,4)>=maxcoord-dx/2
                k = k + 1;
                RemoveVoxel(1,k) = i;
            end
        end
        
        Nodedofs = [EdofElasticity(:,1)';EdofElasticity(:,2)';EdofElasticity(:,3)']';
        EdofElasticity(RemoveVoxel,:) = [];
        
        
        % hack to have 'ndof' as a struct
        ndof =  meshProperties.ndof;
        meshProperties = rmfield(meshProperties,'ndof');
        meshProperties.ndof.diffusion = ndof;
        meshProperties.ndof.elasticity = max(max(EdofElasticity));
        
        
        % Generates sparse matrix structure % K = sparse(sparseVec.i,sparseVec.j,X,ndof,ndof)
        first = 1:floor(meshProperties.nel/4);
        second = first(end)+1:floor(meshProperties.nel/2);
        third = second(end)+1:floor(meshProperties.nel*3/4);
        fourth = third(end)+1:meshProperties.nel;
        
        firstSparseVecElasticity.i = zeros(length(first)*8*3,8*3,'uint32');
        firstSparseVecElasticity.j = zeros(length(first)*8*3,8*3,'uint32');
        firstSparseVecElasticity.first = first;
        
        secondSparseVecElasticity.i = zeros(length(second)*8*3,8*3,'uint32');
        secondSparseVecElasticity.j = zeros(length(second)*8*3,8*3,'uint32');
        secondSparseVecElasticity.second = second;
        
        thirdSparseVecElasticity.i = zeros(length(third)*8*3,8*3,'uint32');
        thirdSparseVecElasticity.j = zeros(length(third)*8*3,8*3,'uint32');
        thirdSparseVecElasticity.third = third;
        
        fourthSparseVecElasticity.i = zeros(length(fourth)*8*3,8*3,'uint32');
        fourthSparseVecElasticity.j = zeros(length(fourth)*8*3,8*3,'uint32');
        fourthSparseVecElasticity.fourth = fourth;
        
        kk = 0;
         for i=first
             % elasticity
             for j=1:8*3
                 kk = kk + 1;
                 firstSparseVecElasticity.i(kk,:) = EdofElasticity(i,j)*ones(1,8*3,'uint32');
                 firstSparseVecElasticity.j(kk,:) = EdofElasticity(i,:);
             end
         end
         
         kk = 0;
         for i=second
             % elasticity
             for j=1:8*3
                 kk = kk + 1;
                 secondSparseVecElasticity.i(kk,:) = EdofElasticity(i,j)*ones(1,8*3,'uint32');
                 secondSparseVecElasticity.j(kk,:) = EdofElasticity(i,:);
             end
         end
         
         kk = 0;
         for i=third
             % elasticity
             for j=1:8*3
                 kk = kk + 1;
                 thirdSparseVecElasticity.i(kk,:) = EdofElasticity(i,j)*ones(1,8*3,'uint32');
                 thirdSparseVecElasticity.j(kk,:) = EdofElasticity(i,:);
             end
         end
         
         kk = 0;
         for i=fourth
             % elasticity
             for j=1:8*3
                 kk = kk + 1;
                 fourthSparseVecElasticity.i(kk,:) = EdofElasticity(i,j)*ones(1,8*3,'uint32');
                 fourthSparseVecElasticity.j(kk,:) = EdofElasticity(i,:);
             end
         end
        
        
        firstSparseVecElasticity.i = reshape(firstSparseVecElasticity.i',length(first)*8*3*8*3,1);
        firstSparseVecElasticity.j = reshape(firstSparseVecElasticity.j',length(first)*8*3*8*3,1);
        secondSparseVecElasticity.i = reshape(secondSparseVecElasticity.i',length(second)*8*3*8*3,1);
        secondSparseVecElasticity.j = reshape(secondSparseVecElasticity.j',length(second)*8*3*8*3,1);
        thirdSparseVecElasticity.i = reshape(thirdSparseVecElasticity.i',length(third)*8*3*8*3,1);
        thirdSparseVecElasticity.j = reshape(thirdSparseVecElasticity.j',length(third)*8*3*8*3,1);
        fourthSparseVecElasticity.i = reshape(fourthSparseVecElasticity.i',length(fourth)*8*3*8*3,1);
        fourthSparseVecElasticity.j = reshape(fourthSparseVecElasticity.j',length(fourth)*8*3*8*3,1);
        
        
        % Determines 2D boundary elements along all faces
        [boundaryDofs] = obj.findBoundaryDofs(NodeCoords,Nodedofs,meshProperties,obj.size);
        
        saveVoxelCoords = [obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'];
        save(saveVoxelCoords,'EdofElasticity','firstSparseVecElasticity','secondSparseVecElasticity','thirdSparseVecElasticity','fourthSparseVecElasticity','meshProperties','boundaryDofs','elementMaterial','Nodedofs','-v7.3','-append');
        
        end
        function getBoundaryElements(obj)
        load([obj.path2Realization,'SVEparameters_',num2str(obj.realizationNumber),'.mat'])
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'Edof')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'meshProperties')
        load([obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'],'NodeCoords')
        
        boundaryElement = zeros(meshProperties.nel,1);
        for iel=1:meshProperties.nel
            elementCenterCoords = NodeCoords(Edof(iel,2),2:end) + [0.5 0.5 0.5]*meshProperties.dx;
            
%             % x
%             if elementCenterCoords(1) < meshProperties.dx || elementCenterCoords(1) > obj.size - meshProperties.dx
%                 boundaryElement(iel) = 1;
%             end
%             
%             % y
%             if elementCenterCoords(2) < meshProperties.dx || elementCenterCoords(2) > obj.size - meshProperties.dx
%                 boundaryElement(iel) = 1;
%             end
            
            % z
            if elementCenterCoords(3) < meshProperties.dx || elementCenterCoords(3) > obj.size - meshProperties.dx
                boundaryElement(iel) = 1;
            end
        end
        saveVoxelCoords = [obj.path2Realization,'TopologyBundle_',num2str(obj.nx),'_',num2str(obj.realizationNumber),'.mat'];
        save(saveVoxelCoords,'boundaryElement','-v7.3','-append');
        end
    end
    methods(Access = private)
        function [interfaceVoxel,Edof] = findITZ(obj,meshProperties,centroid,radius,NodeCoords,Edof)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % findITZ.m
        %
        % Finds line-sphere intersections for all voxels
        % http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
        %
        % Written by Filip Nilenius, 2013-01-10
        % Last edited on: 2013-09-30
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Finds line-sphere intersections for all voxels
        % http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

        normal  = [1 0 0    % line 1
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

        nodePick = 1+[1 2 4 1 5 6 8 5 1 2 3 4];

        k = 0;
        kkk = 0;

        interfaceVoxel.volume.cement = zeros(meshProperties.nel,1);
        interfaceVoxel.volume.ballast = zeros(meshProperties.nel,1);
        interfaceVoxel.area.ITZ = zeros(meshProperties.nel,1);
        interfaceVoxel.surfaceNormal = zeros(meshProperties.nel,3);
        interfaceVoxel.boundaryAggregate = zeros(meshProperties.nel,1);

        for iel=1:meshProperties.nel

            % identify which ballast particle is closest to the current voxel
            voxelCenter = (NodeCoords(Edof(iel,2),2:end) + NodeCoords(Edof(iel,8),2:end))/2;
            AA = bsxfun(@minus,voxelCenter,centroid);
            normA = sqrt(sum(AA.^2,2));
            distanceToSurface = normA - radius;
            [Y closestBallast] = min(distanceToSurface);

            if distanceToSurface(closestBallast) <= radius(closestBallast)+3*meshProperties.dx % only check voxel if reasonably close to ballast surface
                kk = 0;
                elementPointCollector = zeros(1,3);
                for iline = 1:12 % for all 12 line segments of the voxel
                    origin = NodeCoords(Edof(iel,nodePick(iline)),2:4);
                    A = origin - centroid(closestBallast,:);

                    % if any intersection exists
                    apa = (normal(iline,:)*A')^2 - sum(A'.^2) + radius(closestBallast)^2;
                    if apa>=0
                        d(iline).plus  = -normal(iline,:)*A' + sqrt(apa);
                        d(iline).minus = -normal(iline,:)*A' - sqrt(apa);

                        if d(iline).plus >= 0 && d(iline).plus <= meshProperties.dx
                            kk = kk + 1;
                            elementPointCollector(kk,:) = origin + d(iline).plus*normal(iline,:);
                        elseif d(iline).minus >= 0 && d(iline).minus <= meshProperties.dx
                            kk = kk + 1;
                            elementPointCollector(kk,:) = origin + d(iline).minus*normal(iline,:);
                        end
                    end
                end
                
                
                if length(elementPointCollector(:,1))>2 % if more than 2 intersection points exist
                    % determine ballast/cement paste nodes at voxel corners
                    counter.ballast = 0;
                    counter.cement = 0;
                    for i=1:8
                        nodeCoord = NodeCoords(Edof(iel,i+1),2:4);
                        Norm = norm(nodeCoord - centroid(closestBallast,:));

                        if Norm>radius(closestBallast)
                            counter.cement = counter.cement + 1;
                            node.cement(counter.cement,:) = nodeCoord;
                        else
                            counter.ballast = counter.ballast + 1;
                            node.ballast(counter.ballast,:) = nodeCoord;
                        end
                    end

                    % workaround if 3 points end up on 1 voxel surface
                    if counter.cement ~= 8 && counter.ballast ~= 8

                        Edof(iel,1) = 2;

                        node.cement = [node.cement;elementPointCollector];
                        node.ballast = [node.ballast;elementPointCollector];

                        % determines surface normal
                        surfaceNormal = (node.ballast(end,:) - centroid(closestBallast,:))/norm(node.ballast(end,:) - centroid(closestBallast,:));

                        [volume.cement,area.cement]=area3d(node.cement);
                        [volume.ballast,area.ballast]=area3d(node.ballast);

                        area.ITZ = (area.cement + area.ballast - meshProperties.dx^2*6)/2;

                        interfaceVoxel.volume.cement(iel) = volume.cement;
                        interfaceVoxel.volume.ballast(iel) = volume.ballast;
                        interfaceVoxel.area.ITZ(iel) = area.ITZ;
                        interfaceVoxel.surfaceNormal(iel,:) = surfaceNormal;


                        clear node
                    end
                    
                    % check to see if aggregate cuts through SVE boundary
                    if centroid(closestBallast,3) - radius(closestBallast) < 2*meshProperties.dx || centroid(closestBallast,3) + radius(closestBallast) > obj.size - 2*meshProperties.dx
                        interfaceVoxel.boundaryAggregate(iel) = 1;
                    end
                end
            end
        end

        end
        function [ex,ey,ez] = prepareKe(obj,dx)

        p1 = [0 0 0]*dx;
        p2 = [1 0 0]*dx;
        p3 = [1 1 0]*dx;
        p4 = [0 1 0]*dx;
        p5 = [0 0 1]*dx;
        p6 = [1 0 1]*dx;
        p7 = [1 1 1]*dx;
        p8 = [0 1 1]*dx;

        ex = [p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1) p8(1)];
        ey = [p1(2) p2(2) p3(2) p4(2) p5(2) p6(2) p7(2) p8(2)];
        ez = [p1(3) p2(3) p3(3) p4(3) p5(3) p6(3) p7(3) p8(3)];
        end
        function Ce = Ce_3D(obj,ex,ey,ez,ep,Cap)
        % Ke=flw3i8e(ex,ey,ez,ep,D)
        % [Ke,fe]=flw3i8e(ex,ey,ez,ep,D,eq)
        %-------------------------------------------------------------
        % PURPOSE
        %  Compute element stiffness (conductivity)
        %  matrix for 8 node isoparametric field element
        %
        %  INPUT:  ex = [x1 x2 x3 ... x8] 
        %          ey = [y1 y2 y3 ... y8]       element coordinates
        %          ez = [z1 z2 z3 ... z8] 
        %  
        %          ep = [ir]                    Ir: Integration rule
        %
        %          D  = [kxx kxy kxz;
        %                kyx kyy kyz;
        %                kzx kzy kzz]           constitutive matrix
        %
        %          eq                           heat supply per unit 
        %                                       volume  
        %
        %  OUTPUT: Ke :  element 'stiffness' matrix (8 x 8)
        %
        %          fe :  element load vector (8 x 1)
        %-------------------------------------------------------------

        % LAST MODIFIED: K Persson    1995-08-24
        % Copyright (c)  Division of Structural Mechanics and
        %                Department of Solid Mechanics.
        %                Lund Institute of Technology
        %-------------------------------------------------------------
          ir=ep(1);  ngp=ir*ir*ir;
          if nargin==5; eq=0 ; end


          if ir==2
            g1=0.577350269189626; w1=1;
            gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1; w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
            gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1; w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
            gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1; w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
          elseif ir==3
            g1=0.774596669241483; g2=0.;
            w1=0.555555555555555; w2=0.888888888888888;

            I1=[-1; 0; 1;-1; 0; 1;-1; 0; 1]';
            I2=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]';
            gp(:,1)=[I1 I1 I1]'*g1;
            gp(:,1)=[I2 I2 I2]'*g2+gp(:,1);
            I1=abs(I1);
            I2=abs(I2);
            w(:,1)=[I1 I1 I1]'*w1;
            w(:,1)=[I2 I2 I2]'*w2+w(:,1);
            I1=[-1;-1;-1; 0; 0; 0; 1; 1; 1]';
            I2=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]';
            gp(:,2)=[I1 I1 I1]'*g1;
            gp(:,2)=[I2 I2 I2]'*g2+gp(:,2);
            I1=abs(I1);
            I2=abs(I2);
            w(:,2)=[I1 I1 I1]'*w1;
            w(:,2)=[I2 I2 I2]'*w2+w(:,2);
            I1=[-1;-1;-1;-1;-1;-1;-1;-1;-1]';
            I2=[ 0; 0; 0; 0; 0; 0; 0; 0; 0]';
            I3=abs(I1);
            gp(:,3)=[I1 I2 I3]'*g1;
            gp(:,3)=[I2 I3 I2]'*g2+gp(:,3);
            w(:,3)=[I3 I2 I3]'*w1;
            w(:,3)=[I2 I3 I2]'*w2+w(:,3);
          else
            disp('Used number of integration points not implemented');
            return
          end

          wp=w(:,1).*w(:,2).*w(:,3);


          xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);  r2=ngp*3;

          N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;  N(:,5)=(1-xsi).*(1-eta).*(1+zet)/8;
          N(:,2)=(1+xsi).*(1-eta).*(1-zet)/8;  N(:,6)=(1+xsi).*(1-eta).*(1+zet)/8;
          N(:,3)=(1+xsi).*(1+eta).*(1-zet)/8;  N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
          N(:,4)=(1-xsi).*(1+eta).*(1-zet)/8;  N(:,8)=(1-xsi).*(1+eta).*(1+zet)/8;

          dNr(1:3:r2,1)=-(1-eta).*(1-zet);    dNr(1:3:r2,2)= (1-eta).*(1-zet);
          dNr(1:3:r2,3)= (1+eta).*(1-zet);    dNr(1:3:r2,4)=-(1+eta).*(1-zet);
          dNr(1:3:r2,5)=-(1-eta).*(1+zet);    dNr(1:3:r2,6)= (1-eta).*(1+zet);
          dNr(1:3:r2,7)= (1+eta).*(1+zet);    dNr(1:3:r2,8)=-(1+eta).*(1+zet);
          dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet);  dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet);
          dNr(2:3:r2+1,3)= (1+xsi).*(1-zet);  dNr(2:3:r2+1,4)= (1-xsi).*(1-zet);
          dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet);  dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet);
          dNr(2:3:r2+1,7)= (1+xsi).*(1+zet);  dNr(2:3:r2+1,8)= (1-xsi).*(1+zet);
          dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta);  dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta);
          dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta);  dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta);
          dNr(3:3:r2+2,5)= (1-xsi).*(1-eta);  dNr(3:3:r2+2,6)= (1+xsi).*(1-eta);
          dNr(3:3:r2+2,7)= (1+xsi).*(1+eta);  dNr(3:3:r2+2,8)= (1-xsi).*(1+eta);
          dNr=dNr/8.;


          Ce = zeros(8,8);
          JT=dNr*[ex;ey;ez]';

          for i=1:ngp
            indx=[ 3*i-2; 3*i-1; 3*i ];
            detJ=det(JT(indx,:));
            if detJ<10*eps
              disp('Jacobideterminanten lika med noll!')
            end
            Ce = Ce + N(i,:)'*N(i,:)*detJ*wp(i);
          end

          Ce = Ce*Cap;
        end
        function [Ke,fe]=flw3i8e(obj,ex,ey,ez,ep,D,eq)
        % Ke=flw3i8e(ex,ey,ez,ep,D)
        % [Ke,fe]=flw3i8e(ex,ey,ez,ep,D,eq)
        %-------------------------------------------------------------
        % PURPOSE
        %  Compute element stiffness (conductivity)
        %  matrix for 8 node isoparametric field element
        %
        %  INPUT:  ex = [x1 x2 x3 ... x8] 
        %          ey = [y1 y2 y3 ... y8]       element coordinates
        %          ez = [z1 z2 z3 ... z8] 
        %  
        %          ep = [ir]                    Ir: Integration rule
        %
        %          D  = [kxx kxy kxz;
        %                kyx kyy kyz;
        %                kzx kzy kzz]           constitutive matrix
        %
        %          eq                           heat supply per unit 
        %                                       volume  
        %
        %  OUTPUT: Ke :  element 'stiffness' matrix (8 x 8)
        %
        %          fe :  element load vector (8 x 1)
        %-------------------------------------------------------------

        % LAST MODIFIED: K Persson    1995-08-24
        % Copyright (c)  Division of Structural Mechanics and
        %                Department of Solid Mechanics.
        %                Lund Institute of Technology
        %-------------------------------------------------------------
          ir=ep(1);  ngp=ir*ir*ir;
          if nargin==5; eq=0 ; end


          if ir==2
            g1=0.577350269189626; w1=1;
            gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1; w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
            gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1; w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
            gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1; w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
          elseif ir==3
            g1=0.774596669241483; g2=0.;
            w1=0.555555555555555; w2=0.888888888888888;

            I1=[-1; 0; 1;-1; 0; 1;-1; 0; 1]';
            I2=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]';
            gp(:,1)=[I1 I1 I1]'*g1;
            gp(:,1)=[I2 I2 I2]'*g2+gp(:,1)
            I1=abs(I1);
            I2=abs(I2);
            w(:,1)=[I1 I1 I1]'*w1;
            w(:,1)=[I2 I2 I2]'*w2+w(:,1);
            I1=[-1;-1;-1; 0; 0; 0; 1; 1; 1]';
            I2=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]';
            gp(:,2)=[I1 I1 I1]'*g1;
            gp(:,2)=[I2 I2 I2]'*g2+gp(:,2);
            I1=abs(I1);
            I2=abs(I2);
            w(:,2)=[I1 I1 I1]'*w1;
            w(:,2)=[I2 I2 I2]'*w2+w(:,2);
            I1=[-1;-1;-1;-1;-1;-1;-1;-1;-1]';
            I2=[ 0; 0; 0; 0; 0; 0; 0; 0; 0]';
            I3=abs(I1);
            gp(:,3)=[I1 I2 I3]'*g1;
            gp(:,3)=[I2 I3 I2]'*g2+gp(:,3);
            w(:,3)=[I3 I2 I3]'*w1;
            w(:,3)=[I2 I3 I2]'*w2+w(:,3);
          else
            disp('Used number of integration points not implemented');
            return
          end

          wp=w(:,1).*w(:,2).*w(:,3);


          xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);  r2=ngp*3;

          N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;  N(:,5)=(1-xsi).*(1-eta).*(1+zet)/8;
          N(:,2)=(1+xsi).*(1-eta).*(1-zet)/8;  N(:,6)=(1+xsi).*(1-eta).*(1+zet)/8;
          N(:,3)=(1+xsi).*(1+eta).*(1-zet)/8;  N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
          N(:,4)=(1-xsi).*(1+eta).*(1-zet)/8;  N(:,8)=(1-xsi).*(1+eta).*(1+zet)/8;

          dNr(1:3:r2,1)=-(1-eta).*(1-zet);    dNr(1:3:r2,2)= (1-eta).*(1-zet);
          dNr(1:3:r2,3)= (1+eta).*(1-zet);    dNr(1:3:r2,4)=-(1+eta).*(1-zet);
          dNr(1:3:r2,5)=-(1-eta).*(1+zet);    dNr(1:3:r2,6)= (1-eta).*(1+zet);
          dNr(1:3:r2,7)= (1+eta).*(1+zet);    dNr(1:3:r2,8)=-(1+eta).*(1+zet);
          dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet);  dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet);
          dNr(2:3:r2+1,3)= (1+xsi).*(1-zet);  dNr(2:3:r2+1,4)= (1-xsi).*(1-zet);
          dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet);  dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet);
          dNr(2:3:r2+1,7)= (1+xsi).*(1+zet);  dNr(2:3:r2+1,8)= (1-xsi).*(1+zet);
          dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta);  dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta);
          dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta);  dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta);
          dNr(3:3:r2+2,5)= (1-xsi).*(1-eta);  dNr(3:3:r2+2,6)= (1+xsi).*(1-eta);
          dNr(3:3:r2+2,7)= (1+xsi).*(1+eta);  dNr(3:3:r2+2,8)= (1-xsi).*(1+eta);
          dNr=dNr/8.;


          Ke1=zeros(8,8);  fe1=zeros(8,1);
          JT=dNr*[ex;ey;ez]';

          for i=1:ngp
            indx=[ 3*i-2; 3*i-1; 3*i ];
            detJ=det(JT(indx,:));
            if detJ<10*eps
              disp('Jacobideterminanten lika med noll!')
            end
            JTinv=inv(JT(indx,:));
            B=JTinv*dNr(indx,:);
            Ke1=Ke1+B'*D*B*detJ*wp(i);
            fe1=fe1+N(i,:)'*detJ*wp(i);
          end

          Ke=Ke1; fe=fe1*eq;
        end
        function [BEdof] = boundaryEdof(obj,Edof,meshProperties,NodeCoords,Lbox)
            
            BEdof.x.back.Edof  = zeros((meshProperties.nx-1)^2,5);
            BEdof.x.front.Edof = zeros((meshProperties.nx-1)^2,5);
            BEdof.y.back.Edof  = zeros((meshProperties.nx-1)^2,5);
            BEdof.y.front.Edof = zeros((meshProperties.nx-1)^2,5);
            BEdof.z.back.Edof  = zeros((meshProperties.nx-1)^2,5);
            BEdof.z.front.Edof = zeros((meshProperties.nx-1)^2,5);

            BEdof.x.back.counter  = zeros((meshProperties.nx-1)^2,5);
            BEdof.x.front.counter = zeros((meshProperties.nx-1)^2,5);
            BEdof.y.back.counter  = zeros((meshProperties.nx-1)^2,5);
            BEdof.y.front.counter = zeros((meshProperties.nx-1)^2,5);
            BEdof.z.back.counter  = zeros((meshProperties.nx-1)^2,5);
            BEdof.z.front.counter = zeros((meshProperties.nx-1)^2,5);


            for i=1:meshProperties.nel
                % Face BEdof.x.back.Edof
                if NodeCoords(Edof(i,2),2) < meshProperties.dx/2 % tests node 1
                    BEdof.x.back.counter = BEdof.x.back.counter + 1;
                    BEdof.x.back.Edof(BEdof.x.back.counter,1) = Edof(i,2);
                    BEdof.x.back.Edof(BEdof.x.back.counter,2) = Edof(i,6);
                    BEdof.x.back.Edof(BEdof.x.back.counter,3) = Edof(i,9);
                    BEdof.x.back.Edof(BEdof.x.back.counter,4) = Edof(i,5);
                    BEdof.x.back.Edof(BEdof.x.back.counter,5) = Edof(i,1);
                end

                % Face BEdof.x.front.Edof
                if NodeCoords(Edof(i,3),2) > Lbox-meshProperties.dx/2 % tests node 2
                    BEdof.x.front.counter = BEdof.x.front.counter + 1;
                    BEdof.x.front.Edof(BEdof.x.front.counter,1) = Edof(i,7);
                    BEdof.x.front.Edof(BEdof.x.front.counter,2) = Edof(i,3);
                    BEdof.x.front.Edof(BEdof.x.front.counter,3) = Edof(i,4);
                    BEdof.x.front.Edof(BEdof.x.front.counter,4) = Edof(i,8);
                    BEdof.x.front.Edof(BEdof.x.front.counter,5) = Edof(i,1);
                end

                % Face BEdof.y.back.Edof
                if NodeCoords(Edof(i,2),3) < meshProperties.dx/2 % tests node 1
                    BEdof.y.back.counter = BEdof.y.back.counter + 1;
                    BEdof.y.back.Edof(BEdof.y.back.counter,1) = Edof(i,2);
                    BEdof.y.back.Edof(BEdof.y.back.counter,2) = Edof(i,3);
                    BEdof.y.back.Edof(BEdof.y.back.counter,3) = Edof(i,7);
                    BEdof.y.back.Edof(BEdof.y.back.counter,4) = Edof(i,6);
                    BEdof.y.back.Edof(BEdof.y.back.counter,5) = Edof(i,1);
                end

                % Face BEdof.y.front.Edof
                if NodeCoords(Edof(i,5),3) > Lbox-meshProperties.dx/2 % tests node 4
                    BEdof.y.front.counter = BEdof.y.front.counter + 1;
                    BEdof.y.front.Edof(BEdof.y.front.counter,1) = Edof(i,4);
                    BEdof.y.front.Edof(BEdof.y.front.counter,2) = Edof(i,6);
                    BEdof.y.front.Edof(BEdof.y.front.counter,3) = Edof(i,9);
                    BEdof.y.front.Edof(BEdof.y.front.counter,4) = Edof(i,8);
                    BEdof.y.front.Edof(BEdof.y.front.counter,5) = Edof(i,1);
                end

                % Face BEdof.z.back.Edof
                if NodeCoords(Edof(i,2),4) < meshProperties.dx/2 % tests node 1
                    BEdof.z.back.counter = BEdof.z.back.counter + 1;
                    BEdof.z.back.Edof(BEdof.z.back.counter,1) = Edof(i,3);
                    BEdof.z.back.Edof(BEdof.z.back.counter,2) = Edof(i,2);
                    BEdof.z.back.Edof(BEdof.z.back.counter,3) = Edof(i,5);
                    BEdof.z.back.Edof(BEdof.z.back.counter,4) = Edof(i,4);
                    BEdof.z.back.Edof(BEdof.z.back.counter,5) = Edof(i,1);
                end
                
                % Face BEdof.z.front.Edof
                if NodeCoords(Edof(i,6),4) > Lbox-meshProperties.dx/2 % tests node 5
                    BEdof.z.front.counter = BEdof.z.front.counter + 1;
                    BEdof.z.front.Edof(BEdof.z.front.counter,1) = Edof(i,6);
                    BEdof.z.front.Edof(BEdof.z.front.counter,2) = Edof(i,7);
                    BEdof.z.front.Edof(BEdof.z.front.counter,3) = Edof(i,8);
                    BEdof.z.front.Edof(BEdof.z.front.counter,4) = Edof(i,9);
                    BEdof.z.front.Edof(BEdof.z.front.counter,5) = Edof(i,1);
                end
            end
        end
        function [boundaryDofs] = findBoundaryDofs(obj,NodeCoords,Nodedofs,meshProperties,Lbox)
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here

        % x
        boundaryDofs.x.back.counter = 0;
        boundaryDofs.x.back.nodes = zeros(meshProperties.nx^2,3);
        boundaryDofs.x.front.counter = 0;
        boundaryDofs.x.front.nodes = zeros(meshProperties.nx^2,3);

        % y
        boundaryDofs.y.back.counter = 0;
        boundaryDofs.y.back.nodes = zeros(meshProperties.nx^2,3);
        boundaryDofs.y.front.counter = 0;
        boundaryDofs.y.front.nodes = zeros(meshProperties.nx^2,3);

        % z
        boundaryDofs.z.back.counter = 0;
        boundaryDofs.z.back.nodes = zeros(meshProperties.nx^2,3);
        boundaryDofs.z.front.counter = 0;
        boundaryDofs.z.front.nodes = zeros(meshProperties.nx^2,3);
        
        % line 1
        counter = 0;
        boundaryDofs.line1 = zeros(meshProperties.nx,1);


        for inode = 1:length(NodeCoords)
            % x back
            if NodeCoords(inode,2) < meshProperties.dx/2
                boundaryDofs.x.back.counter = boundaryDofs.x.back.counter + 1;
                boundaryDofs.x.back.nodes(boundaryDofs.x.back.counter,:) = Nodedofs(inode,:);
            end

            % x front
            if NodeCoords(inode,2) > Lbox - meshProperties.dx/2
                boundaryDofs.x.front.counter = boundaryDofs.x.front.counter + 1;
                boundaryDofs.x.front.nodes(boundaryDofs.x.front.counter,:) = Nodedofs(inode,:);
            end

            % y back
            if NodeCoords(inode,3) < meshProperties.dx/2
                boundaryDofs.y.back.counter = boundaryDofs.y.back.counter + 1;
                boundaryDofs.y.back.nodes(boundaryDofs.y.back.counter,:) = Nodedofs(inode,:);
            end

            % y front
            if NodeCoords(inode,3) > Lbox - meshProperties.dx/2
                boundaryDofs.y.front.counter = boundaryDofs.y.front.counter + 1;
                boundaryDofs.y.front.nodes(boundaryDofs.y.front.counter,:) = Nodedofs(inode,:);
            end

            % z back
            if NodeCoords(inode,4) < meshProperties.dx/2
                boundaryDofs.z.back.counter = boundaryDofs.z.back.counter + 1;
                boundaryDofs.z.back.nodes(boundaryDofs.z.back.counter,:) = Nodedofs(inode,:);
            end

            % z front
            if NodeCoords(inode,4) > Lbox - meshProperties.dx/2
                boundaryDofs.z.front.counter = boundaryDofs.z.front.counter + 1;
                boundaryDofs.z.front.nodes(boundaryDofs.z.front.counter,:) = Nodedofs(inode,:);
            end
            
            % line 1
            if NodeCoords(inode,2) < meshProperties.dx/2 && NodeCoords(inode,3) < meshProperties.dx/2
                counter = counter + 1;
                boundaryDofs.line1(counter) = Nodedofs(inode,2);
            end
        end
        end
        function makeMovie(obj,v,gravelSet,cube)
            % make movie
            coords = reshape([gravelSet.coordinates],3,length(gravelSet))';
            hold on
            bubbleplot3(coords(:,1),coords(:,2),coords(:,3),[gravelSet.radius]');
            plot3([0 obj.Lbox],[0 0],[0 0],'color',[0.3 0.3 0.3])
            plot3([0 0],[0 obj.Lbox],[0 0],'color',[0.3 0.3 0.3])
            plot3([0 0],[0 0],[0 obj.Lbox],'color',[0.3 0.3 0.3])
            
            plot3([obj.Lbox obj.Lbox],[obj.Lbox obj.Lbox],[0 obj.Lbox],'color',[0.3 0.3 0.3])
            plot3([obj.Lbox 0],[obj.Lbox obj.Lbox],[0 0],'color',[0.3 0.3 0.3])
            plot3([obj.Lbox obj.Lbox],[obj.Lbox 0],[0 0],'color',[0.3 0.3 0.3])
            
            plot3([obj.Lbox obj.Lbox],[0 0],[obj.Lbox 0],'color',[0.3 0.3 0.3])
            plot3([obj.Lbox 0],[0 0],[obj.Lbox obj.Lbox],'color',[0.3 0.3 0.3])
            plot3([obj.Lbox obj.Lbox],[0 obj.Lbox],[obj.Lbox obj.Lbox],'color',[0.3 0.3 0.3])
            
            plot3([0 0],[obj.Lbox obj.Lbox],[obj.Lbox 0],'color',[0.3 0.3 0.3])
            plot3([0 0],[obj.Lbox 0],[obj.Lbox obj.Lbox],'color',[0.3 0.3 0.3])
            plot3([0 obj.Lbox],[obj.Lbox obj.Lbox],[obj.Lbox obj.Lbox],'color',[0.3 0.3 0.3])
            
            
            plot3([0 cube.size],[0 0],[0 0],'color','black','LineWidth',2)
            plot3([0 0],[0 cube.size],[0 0],'color','black','LineWidth',2)
            plot3([0 0],[0 0],[0 cube.size],'color','black','LineWidth',2)
            
            plot3([cube.size cube.size],[cube.size cube.size],[0 cube.size],'color','black','LineWidth',2)
            plot3([cube.size 0],[cube.size cube.size],[0 0],'color','black','LineWidth',2)
            plot3([cube.size cube.size],[cube.size 0],[0 0],'color','black','LineWidth',2)
            
            plot3([cube.size cube.size],[0 0],[cube.size 0],'color','black','LineWidth',2)
            plot3([cube.size 0],[0 0],[cube.size cube.size],'color','black','LineWidth',2)
            plot3([cube.size cube.size],[0 cube.size],[cube.size cube.size],'color','black','LineWidth',2)
            
            plot3([0 0],[cube.size cube.size],[cube.size 0],'color','black','LineWidth',2)
            plot3([0 0],[cube.size 0],[cube.size cube.size],'color','black','LineWidth',2)
            plot3([0 cube.size],[cube.size cube.size],[cube.size cube.size],'color','black','LineWidth',2)
            axis([0 obj.Lbox 0 obj.Lbox 0 obj.Lbox])
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            set(gca,'ztick',[])
            set(gca,'zticklabel',[])
            set(gcf,'position',[0,0,720,720])
            
            axis square
            hold off
            frame = getframe(gcf);
            writeVideo(v,frame);
            clf
        end
    end
end