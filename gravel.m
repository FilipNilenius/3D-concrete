classdef gravel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius
        coordinates
        velocity
    end
    properties (Dependent)
        mass
    end
    
    methods
        function a = get.mass(obj)
            a = 4/3*pi*obj.radius^3;
        end
    end
    
end

