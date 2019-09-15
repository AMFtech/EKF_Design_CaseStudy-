
clc
clear all;
close all;

Ts = 0.1;                   % sampling time, dt = 0.1 sec
load object_trajectory.mat
RealPos = x_traj_pos;
[m, n] = size(RealPos);
X_RealPos = RealPos(1,:);
Y_RealPos = RealPos(2,:);

% Observation Covariance Matrix
Sig_X   = 5;   % in radian
Sig_Y   = 5;            % in meter

q = 8 *10^-2;
Ts = 0.1;           % sampling time, dt = 0.1 sec

%% Generate random trajectories

% xPosZ  = X_RealPos + randn(1, n) * Sig_X;
% yPosZ  = Y_RealPos + randn(1, n) * Sig_Y;

xPosZ  = X_RealPos + randn(1, n) * Sig_X;
yPosZ  = Y_RealPos + randn(1, n) * Sig_Y;

figure;
plot(xPosZ, yPosZ, '.', 'color','red', 'LineWidth', 1 );
hold on;
plot(X_RealPos, Y_RealPos, 'b', 'LineWidth', 3 );

%% EKF operation ----------------------------------------------------------

% Matrices Initialization -----------------------------

% State Transition Matrix
F = [ 1, 0, Ts, 0;
      0, 1, 0   Ts; 
      0, 0, 1,  0
      0, 0, 0,  1]

% Initial State of State Covariance Matrix
Pk_1 = [ 0, 0, 0, 0;
         0, 0, 0, 0;
         0, 0, 0, 0;
         0, 0, 0, 0;]

% Covariance Matrix of process noise... (Velocity Noise)
Q = q * [ Ts^3/3, 0,      Ts^2/2, 0;
          0,      Ts^3/3, 0,      Ts^2/2; 
          Ts^2/2, 0       Ts,     0;
          0,      Ts^2/2, 0       Ts]
      
% Q = q * [ Ts^4/4,  Ts^3/2,  Ts^4/4,  Ts^3/2;
%           Ts^3/2,  Ts^2/1,  Ts^3/2,  Ts^2/1; 
%           Ts^4/4,  Ts^3/2,  Ts^4/4,  Ts^3/2;
%           Ts^3/2,  Ts^2/1,  Ts^3/2,  Ts^2/1 ]
  
% Observation Covariance Matrix
R = [ Sig_X^2, 0;  
      0,         Sig_Y^2]
  
Hk = [1, 0, 0, 0 ;
     0, 1, 0, 0]

% Initial State Vector
Xk_1 = [50; 50; 0; 0]

% Stored State Estimate, outcome of EKF
X_EKF = zeros(4, n);

% Recursive Kalman Operation-----------------------------------------------
for i = 1 : n
    % Prediction phase --------------------------------
    Xp = F*Xk_1;            % Initial state estimate
    Pp = F*Pk_1*F' + Q;     % Initial Process Covariance Matrix 
 
    % Measurement update ------------------------------
    K = (Pp*Hk')/(Hk*Pp*Hk' + R);      % Kalman Gain
    dZk = [xPosZ(1,i); yPosZ(1,i)] - Hk*Xp;
 
    Xk = Xp + K * dZk;
    Pk = Pp - K * Hk * Pp;
    X_EKF(:,i) = Xk;
    
    Xk_1 = Xk;
    Pk_1 = Pk;
    
end

%% Plot figures
figure;
plot( RealPos(1,:), RealPos(2,:), 'color', 'blue', 'LineWidth', 3 );
hold on;
plot( X_EKF(1,:), X_EKF(2,:), ':', 'color', 'red', 'MarkerSize', 2, 'LineWidth', 1 );
figure,
plot( (X_EKF(3,:).^2 + X_EKF(4,:).^2).^0.5, ':', 'color', 'red', 'MarkerSize', 2, 'LineWidth', 1 );
Hk*Pk*Hk'
R


