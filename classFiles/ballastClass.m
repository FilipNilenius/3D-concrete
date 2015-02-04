classdef ballastClass < materialClass
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = ballastClass()
            obj.YoungsModulus = 70e5; % N/(cm^2)
            obj.PoissonsRatio = 0.2;
            obj.diffusionCoefficient = 0;
            
            obj = constructMaterial(obj);
        end
        
    end
    
end
