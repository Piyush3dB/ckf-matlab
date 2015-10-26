classdef CKF
    
    properties
        
        Q;
        R;
        
        nx;
        
        xkk1;
        Pkk1;
        
        xkk;
        Pkk;
        
        
        
    end
    
    methods
        
        % Constructor
        function obj = CKF(Q, R, nx)
            obj.Q = Q;
            obj.R = R;
            obj.nx = nx;
        end
        
        
        % Propagation
        function obj = Predict(obj, xkk, Pkk)
            
            obj.xkk1 = xkk;
            
            obj.Pkk1 = Pkk + obj.Q;
        end
        
        % Correction
        function obj = Update(obj, z)
            
            R = obj.R;
            xkk1 = obj.xkk1;
            Pkk1 = obj.Pkk1;
            
            %%%========================================================================
            %%% Genrate a set of Cubature Points
            %%%========================================================================
            
            nx = obj.nx; % state vector dimension
            
            nPts = 2*nx;   % No. of Cubature Points
            
            CPtArray = sqrt(nPts/2)*[eye(nx) -eye(nx)];
            
            %%%========================================================================
            %%% find the square-root of Pkk1 using Singular Value Decomposition (SVD)
            %%%========================================================================
            
            Pkk1 = 0.5*(Pkk1+Pkk1'); % make Pkk1 to be symmetric (Property of a covariance matrix !)
            
            [U, S, V] =  svd(Pkk1);
            
            Skk1 = 0.5*(U+V)*sqrt(S);
            
            %%%========================================================================
            
            Xi =  repmat(xkk1,1,nPts) + Skk1*CPtArray;
            
            Zi = MstEq(Xi);
            
            zkk1 = sum(Zi,2)/nPts;      % predicted Measurement
            
            X = (Xi-repmat(xkk1,1,nPts))/sqrt(nPts);
            
            Z = (Zi-repmat(zkk1,1,nPts))/sqrt(nPts);
            
            Pzz = Z*Z'+ R; % Innovations Covariance
            
            Pxz = X*Z'; % cross-covariance
            
            G = Pxz*pinv(Pzz);         % Cubature Kalman Gain
            
            xkk = xkk1 + G*(z - zkk1);
            
            Pkk = Pkk1 - G*Pzz*G';
            
            % Save objects
            obj.xkk = xkk;
            obj.Pkk = Pkk;
            
        end
    end
    
    
    
    methods(Static)
        
        function z = MstEq(x)
            
            z = [ 0.8*cos(x(1,:)) - 0.2*cos(x(1,:) + x(2,:)) ;
                0.8*sin(x(1,:)) - 0.2*sin(x(1,:) + x(2,:)) ];
        end
    end
    
end