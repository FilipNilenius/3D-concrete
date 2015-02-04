classdef materialClass
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties
        stiffnessTensor
    end
    
    properties (SetAccess = protected)
        YoungsModulus
        PoissonsRatio
        diffusionCoefficient
        crackDiffusivity = 500000; %500000
    end
    
    properties (SetAccess = private)
        lambda
        shearModulus
    end
    
    methods (Access = protected)
        function obj = constructMaterial(obj)
            obj = getlambda(obj);
            obj = getshearModulus(obj);
            obj = getstiffnessTensor(obj);
            obj = getstiffnessVoigt(obj);
        end
    end
    methods (Access = private)
        function obj = getlambda(obj)
            obj.lambda = obj.YoungsModulus*obj.PoissonsRatio/((1 + obj.PoissonsRatio)*(1-2*obj.PoissonsRatio));
        end
        function obj = getshearModulus(obj)
            obj.shearModulus = obj.YoungsModulus/(2*(1 + obj.PoissonsRatio));
        end
        function obj = getstiffnessTensor(obj)
            I.full = zeros(3,3,3,3);
            I.sym.first = zeros(3,3,3,3);
            I.sym.second = zeros(3,3,3,3);
            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                            % full
                            if i==j && k==l
                                I.full(i,j,k,l) = 1;
                            end
                        end
                    end
                end
            end
            
            for i=1:3
                for k=1:3
                    for j=1:3
                        for l=1:3
                            % sym first term
                            if i==k && j==l
                                I.sym.first(i,j,k,l) = 1;
                            end
                        end
                    end
                end
            end
            
            for i=1:3
                for l=1:3
                    for j=1:3
                        for k=1:3
                            % sym first term
                            if i==l && j==k
                                I.sym.second(i,j,k,l) = 1;
                            end
                        end
                    end
                end
            end
            
            obj.stiffnessTensor.full = obj.lambda*I.full + obj.shearModulus*(I.sym.first + I.sym.second);
        end
        function obj = getstiffnessVoigt(obj)
            obj.stiffnessTensor.Voigt = VoigtMap(obj.stiffnessTensor.full,'tensor2voigt');
        end
    end
    
end

