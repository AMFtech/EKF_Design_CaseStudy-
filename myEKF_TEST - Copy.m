
clc
clear all;
close all;

load object_trajectory.mat
RealPos = x_traj_pos;
[m, n] = size(RealPos);
X_RealPos = RealPos(1,:);
Y_RealPos = RealPos(2,:);

% Observation Covariance Matrix
Sig_theta = 3 * (pi/180);   % in radian
Sig_rho   = 0.1;            % in meter

q = 5;
Ts = 0.1;           % sampling time, dt = 0.1 sec

%% Generate random trajectories
[THETA, RHO] = cart2pol(X_RealPos, Y_RealPos);
THETAz  = THETA + randn(1, n) * Sig_theta;
RHOz    = RHO   + randn(1, n) * Sig_rho;
[xPosZ, yPosZ] = pol2cart(THETAz, RHOz);
Z_noise = [X_RealPos ; Y_RealPos] - [xPosZ ; yPosZ];

figure;
plot(xPosZ, yPosZ, '--.', 'color','red', 'LineWidth', 1 );
hold on;
plot(X_RealPos, Y_RealPos, 'b', 'LineWidth', 3 );
    Str = ['Red: Sonar Data, Blue: Real Trajectory, q = ' num2str(q)];
    title(Str, 'FontSize', 12, ...
        'FontWeight','bold');
    xlabel('X position in Cartesian Coordinate', 'FontSize', 10, ...
        'FontWeight','bold');
    ylabel('Y position in Cartesian Coordinate', 'FontSize', 10, ...
        'FontWeight','bold');  
    
% figure;
% plot(1:1:n, THETA*180/pi); hold on; plot(1:1:n, THETAz*180/pi,'--r');
% figure;
% plot(RHO); hold on; plot(RHOz,'--r');


% figure;
% plot(THETA,'--b', 'LineWidth', 2); 
% hold on; plot(THETAz,'.', 'LineWidth', 2); 
% 
% figure;
% plot(Z_noise(1,:)./X_RealPos,'.', 'LineWidth', 2); 
% figure;
% plot(Z_noise(2,:)./Y_RealPos,'.', 'LineWidth', 2); 

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
  
% Observation Covariance Matrix
R = [ Sig_rho^2, 0;  
      0,         Sig_theta^2]

% Initial State Vector
Xk_1 = [50; 50; 0; 0]

% Stored State Estimate, outcome of EKF
X_EKF = zeros(4, n);

% Recursive Kalman Operation-----------------------------------------------
for i = 1 : 500
    % Prediction phase --------------------------------
    Xp = F*Xk_1;            % Initial state estimate
    Pp = F*Pk_1*F' + Q;     % Initial Process Covariance Matrix 
 
    % Measurement update ------------------------------
    [THETAp, RHOp] = cart2pol(Xp(1), Xp(2));
    Hk = [cos(THETAp),       sin(THETAp),      0, 0;
          -sin(THETAp)/RHOp, cos(THETAp)/RHOp, 0, 0];
    K = (Pp*Hk')/(Hk*Pp*Hk' + R);      % Kalman Gain
    dZk = [RHOz(1,i); THETAz(1,i)] - [RHOp ; THETAp] ;
 
    Xk = Xp + K * dZk;
    Pk = Pp - K * Hk * Pp;
    X_EKF(:,i) = Xk;
    
    Xk_1 = Xk;
    Pk_1 = Pk;
    
    K_factor(1,i) = norm( Hk*Pp*Hk' / R );
    
end

%% Plot figures

figure;
plot(K_factor,'--.');

figure;
plot( RealPos(1,:), RealPos(2,:), 'color', 'blue', 'LineWidth', 3 );
hold on;
plot( X_EKF(1,:), X_EKF(2,:), ':', 'color', 'red', 'MarkerSize', 2, 'LineWidth', 1 );
figure,
plot( (X_EKF(3,:).^2 + X_EKF(4,:).^2).^0.5, ':', 'color', 'red', 'MarkerSize', 2, 'LineWidth', 1 );

