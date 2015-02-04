clear all
clc

addpath '\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\solvers';
addpath '\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\postprocessor';
addpath '\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\misc';

LboxID = 6;

for aggFrac = [0.5 0.6 0.65]
    if aggFrac == 0.1
        nx = 180;
    elseif aggFrac == 0.2
        nx = 180;
    elseif aggFrac == 0.3
        nx = 180;
    elseif aggFrac == 0.4
        nx = 200;
    elseif aggFrac == 0.5
        nx = 220;
    elseif aggFrac == 0.6
        nx = 220;
    elseif aggFrac == 0.65
        nx = 220;
    end
    
    for iRealization = 1:20
        path2Realization = ['\\sol.ita.chalmers.se\groups\bom-kt\vc_sem\Betongbyggnad_KG\Projekt\1047Klorid_transport\ber&mod\MATLAB\3D\SVEs\Lbox_',num2str(LboxID),'\aggFrac_',num2str(100*aggFrac),'\Realization_' num2str(iRealization)];
        addpath(path2Realization);
        
        load(['TopologyBundle_',num2str(nx),'_',num2str(iRealization),'.mat'],'Edof')
        load(['TopologyBundle_',num2str(nx),'_',num2str(iRealization),'.mat'],'interfaceVoxel')
        load(['TopologyBundle_',num2str(nx),'_',num2str(iRealization),'.mat'],'meshProperties')
        
        DTaylor = 0;
        for iel=1:meshProperties.nel
            [D] = constitutiveModel(Edof(iel,1),iel,interfaceVoxel);
            DTaylor = DTaylor + D.voxel(1,1);
        end
        DTaylor = DTaylor/meshProperties.nel;
        dlmwrite(['DTaylor_ITZ=000_',num2str(iRealization),'.txt'], DTaylor, 'delimiter', '\t');
        movefile(['DTaylor_ITZ=000_',num2str(iRealization),'.txt'],path2Realization);
        
        % Clear SVE dependencies
        clear Edof
        clear interfaceVoxel
        clear meshProperties
        clear path2Realization
    end
end

