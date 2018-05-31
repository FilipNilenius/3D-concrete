function Ce = mass3i8e(ex,ey,ez,Cap)
            % this function is entirely based on flw3i8e.m, belonging to the
            % CALFEM library, with a few modifications.
            %
            %
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
            ir=3;  ngp=ir*ir*ir;
            if nargin==5; eq=0 ; end
            
            gp = zeros(27,3);
            w = zeros(27,3);
            N = zeros(27,8);
            dNr = zeros(81,8);

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
            Ce = Ce + N(i,:)'*N(i,:)*detJ*wp(i);
            end

            Ce = Ce*Cap;
        end