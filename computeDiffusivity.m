function computeDiffusivity()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computeDiffusivity.m
%
% 'computeDiffusivity.m' wraps around LinStatSolver.m and
% StatPostProcessor.m to compute components of the homogenized
% diffusivity tensor. The implementation is parallelized, i.e. the number 
% of avaiable cores are the number of SVE problems that can be solved
% in paralell.
%
% Written by Filip Nilenius, 2012-06-18
% Last edited on: 2013-03-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

location = 1;

if location == 1 % bom server
	addpath '\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\solvers';
	addpath '\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\postprocessor';
	addpath '\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\misc';
elseif location == 0 % cluster
	addpath '/c3se/users/v03nifi/Glenn/matlab/3D/solvers'
	addpath '/c3se/users/v03nifi/Glenn/matlab/3D/postprocessor'
	addpath '/c3se/users/v03nifi/Glenn/matlab/3D/misc'
end

for LboxID = [2]
    
    % Specify Lbox dependencies
    if LboxID == 2
        nx=70;
        realizationID = [1];
        pauseMinutes = 0;
%         matlabpool open local 10
    elseif LboxID == 4
        nx=120;
        realizationID = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40];
        pauseMinutes = 0;
		nx = 220;
        matlabpool open local 10
    elseif LboxID == 6
		nx = 220;
        realizationID = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    elseif LboxID == 8
        nx=220;
        realizationID = [1 2 3];
        pauseMinutes = 1;
        matlabpool open local 3
    elseif LboxID == 10
        nx=260;
        realizationID = [5 6];
        pauseMinutes = 1;
        matlabpool open local 2
    end
    
%     for aggFrac = [0.60 0.65]
% 		matlabpool open local 10
		for iRealization=1:length(realizationID)
			pause(iRealization*45) % to get the workers out of sync.
			computeDiffusivityFunc(realizationID(iRealization),LboxID,nx,location);
		end
% 		matlabpool close
% 	end
	%rmdir('/c3se/users/v03nifi/Glenn/.matlab','s');
end
end


function computeDiffusivityFunc(iRealization,LboxID,nx,location)

if location == 1 % bom
    path2Realization = ['\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\SVEs\Lbox_',num2str(LboxID),'\Realization_' num2str(iRealization)];
elseif location == 0 % cluster
    path2Realization = ['/c3se/users/v03nifi/Glenn/matlab/3D/SVEs/Lbox_',num2str(LboxID),'/aggFrac_',num2str(100*aggFrac),'/Realization_' num2str(iRealization)]; % cluster
end


homoDiff = zeros(3,3);

% sliceVector = 2:nx-2; % 2D
sliceVector = 1; % 3D

H.bar = 1;
for i=1 % impose macro gradient
    if i==1
        H.grad = [-1  0  0];
    elseif i==2
        H.grad = [ 0 -1  0];
    else
        H.grad = [ 0  0 -1];
    end
    a = LinStatSolver(iRealization,H,path2Realization,nx,sliceVector);
    homoDiff(:,i) = StatPostProcessor(a,iRealization,path2Realization,nx,sliceVector);
end

% dlmwrite(['homoDiffTensor_ITZ=0_',num2str(iRealization),'.txt'], homoDiff, 'delimiter', '\t');
% movefile(['homoDiffTensor_ITZ=0_',num2str(iRealization),'.txt'],path2Realization);
end