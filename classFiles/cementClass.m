classdef cementClass < materialClass
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    properties (SetAccess = private)
        Ec = 35e5;
        crackingStress = 3.5e2; % N/(cm^2)
        crackingStrain = 1e-4;
    end
    
    methods
        function obj = cementClass()
            obj.YoungsModulus = 35e5; % N/(cm^2)
            obj.PoissonsRatio = 0.2;
            obj.diffusionCoefficient = 1;
            
            obj = constructMaterial(obj);
        end
    end
    
end

