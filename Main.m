%%%=======================================================================
%%% This matlab code implements the CKF
%%%=======================================================================

clear;
clc;
close all;

nx = 2; % state vector dimension
nState = nx;

nExpt = 10;             % No. of Experiments/Trials

N = 630;                % No. of Time steps

%%% Process Noise Covariance Q
sigma1 = 1e-2;
sigma2 = 1e-1;
Q = diag([sigma1^2  sigma2^2 ]);

%%% Measurement Noise Covariance R
R = 0.005*eye(nState);

% Initialise metric store
MSE    = zeros(nState, N);
estMSE = zeros(nState, N);

xestArray = zeros(nState, N);

% Initial covariance prior
Skk = diag([0.9  pi/6]);
Pkk = Skk*Skk';

% Generate Testcase
[xArray,zArray] = GenerateScenario(Q);

%% Iterate
for expt = 1:nExpt
    
    % Initial mean prior
    xkk = [0.3+0.9*rand; pi/2+pi*rand];

    fprintf('MC Run in Process = %d\n',expt);
    
    % Iterate for evert measurement
    for k = 1:N
        
        % Propagate estimate and covariance
        [xkk1,Pkk1] = Predict(xkk, Pkk, Q);
        
        % Update estimate and covariance
        z = zArray(:,k); % measurement
        [xkk,Pkk] = Update(xkk1, Pkk1, z, R);
        
        % Save state estimate
        xestArray(:,k) = xkk;
        
        % True value
        xTrue = xArray(:, k);
        
        % Calculate metrics
        MSE(:,k)= MSE(:,k) + sum((xTrue - xkk).^2);
        
        estMSE(:,k) = estMSE(:,k)+ trace(Pkk);
        
    end;    % time-step
    
end;    % expts

MSE    = MSE/(2*nExpt);
estMSE = estMSE/(2*nExpt);

RMSE    = MSE.^(0.5);
estRMSE = estMSE.^(0.5);

%%%========================================================================
%%% Plotting
%%%========================================================================

figure;
subplot(2,1,1);
plot(xArray(1,:),'k');
hold on;
plot(xestArray(1,:),'r:');
ylabel('y(m)','fontsize',16);
legend('Actual','CKF',2);
hold off;

subplot(2,1,2);
plot(xArray(2,:),'k');
hold on;
plot(xestArray(2,:),'r:');
xlabel('Time, k','fontsize',12);
ylabel('y(m)','fontsize',16);
legend('Actual','CKF',2);
hold off;

figure;
x = 1:N;
semilogy(x,RMSE(1,:),'r');
hold on;
semilogy(x,estRMSE(1,:),'r:');
legend('RMSE','filter est. RMSE');
ylabel('RMSE','fontsize',16);
xlabel('Time, k','fontsize',12);


