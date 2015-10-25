%%%=======================================================================
%%% This matlab code implements the CKF 
%%% Author I. Arasaratnam (haran@ieee.org)
%%% Date. May 19, 2012
%%%=======================================================================

clear;
clc; 
close all;

%global Q R;

nExpt = 10;             % No. of Experiments/Trials          

N = 630;                % No. of Time steps            

%%% Process Noise Covariance Q
sigma_processNoise1 = 1e-2;              
sigma_processNoise2 = 1e-1;              
Q = diag([sigma_processNoise1^2  sigma_processNoise2^2 ]); 

%%% Measurement Noise Covariance R
R = 0.005*eye(2);

MSE = zeros(2,N);
estMSE = zeros(2,N);

xestArray = zeros(2,N); 

% Generate Testcase
[xArray,zArray] = GenerateScenario(Q);

for expt = 1:nExpt
    
    %%%====================================================================
    %%% Initialization
    %%%====================================================================
     
    xkk = [0.3+0.9*rand; pi/2+pi*rand];
    
    Skk = diag([0.9  pi/6]);
    Pkk = Skk*Skk';
    
    fprintf('MC Run in Process = %d\n',expt);
    
    for k = 1:N
        
        % Propagate estimate and covariance
        [xkk1,Pkk1] = Predict(xkk,Pkk,Q);
        
        % Update estimate and covariance
        [xkk,Pkk] = Update(xkk1,Pkk1,zArray(:,k),R);
        
        xestArray(:,k) = xkk;
        
        x = xArray(:,k);
        
        MSE(:,k)= MSE(:,k) + sum((x - xkk).^2);
        
        estMSE(:,k) = estMSE(:,k)+ trace(Pkk);
        
    end;    % time-step
    
end;    % expts

MSE = MSE/(2*nExpt);
estMSE = estMSE/(2*nExpt);

RMSE = MSE.^(0.5);
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


