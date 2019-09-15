
clc
clear all;
close all;

load object_trajectory.mat
RealPos = x_traj_pos;
[m, n] = size(RealPos);
X_RealPos = RealPos(1,:);
Y_RealPos = RealPos(2,:);

% Observation Covariance Matrix
global Sig_rho;     % RHO variance in meter
global Sig_theta;   % Theta variance in rad

Sig_theta = 3 * (pi/180);   % in radian
Sig_rho   = 0.1;            % in meter

global q;           % Process Covariance parameter
q = 0.05;
Ts = 0.1;           % sampling time, dt = 0.1 sec

%% Generate a single random trajectory "Zk"

[THETA, RHO] = cart2pol(X_RealPos, Y_RealPos);
THETAz  = THETA + randn(1, n) * Sig_theta;
RHOz    = RHO   + randn(1, n) * Sig_rho;
[xPosZ, yPosZ] = pol2cart(THETAz, RHOz);
Zk_Sonar = [RHOz ; THETAz];

% figure;
% plot(xPosZ, yPosZ, '--.', 'color','green', 'LineWidth', 1 );
% hold on;
% plot(X_RealPos, Y_RealPos, 'b', 'LineWidth', 3 );
%     Str = ['Red: Sonar Data, Blue: Real Trajectory, q = ' num2str(q)];
%     title(Str, 'FontSize', 12, ...
%         'FontWeight','bold');
%     xlabel('X position in Cartesian Coordinate', 'FontSize', 10, ...
%         'FontWeight','bold');
%     ylabel('Y position in Cartesian Coordinate', 'FontSize', 10, ...
%         'FontWeight','bold');  

%% EKF operation ----------------------------------------------------------

[X_EKF, K_Den_Factor, Kk_norm] = myEKF(Zk_Sonar);

%% Plot figures
T = 0:Ts:Ts*(n-1);

figure;
plot(T, xPosZ, 'o', 'color', 'green', 'MarkerSize', 2 );
hold on;
plot(T, yPosZ, 'o', 'color', 'blue', 'MarkerSize', 2 );
hold on;
    legend('The x component', 'The y component', 'Location', 'Best', 'FontSize' , '12', 'FontWeight', 'bold');
    Str = ['The Cartesian components (x & y) of each sensory oservation'];
    title(Str, 'FontSize', 12, 'FontWeight', 'bold');    
    xlabel('Time (sec)', 'FontSize', 10, ...
        'FontWeight','bold');
    ylabel('Cartesin component values of x and y', 'FontSize', 10, ...
        'FontWeight','bold'); 

    % plot(T, K_Den_Factor,'--r', 'LineWidth', 2);
% hold on;
% % plot(T, Kk_norm, 'color', 'blue');
%     Str = ['The Ratio of ||I-K*H)|| / ||K||, q = ', ...
%         num2str(q)];
%     title(Str, 'FontSize', 12, 'FontWeight', 'bold');    
%     xlabel('Time (sec)', 'FontSize', 10, 'FontWeight','bold');
%     ylabel('The Ratio', 'FontSize', 10, 'FontWeight','bold');    
%     grid on;   

figure;
plot(xPosZ, yPosZ, 'o', 'color', 'green', 'MarkerSize', 2 );
hold on;
plot( RealPos(1,:), RealPos(2,:), 'color', 'blue', 'LineWidth', 2 );
hold on;
plot( X_EKF(1,:), X_EKF(2,:), ':', 'color', 'red', 'LineWidth', 3 );
    legend('Sonar Data', 'Real Trajectory', 'EKF model', 'Location', 'Best');
    Str = ['EKF Trajectory compared with the Real Trajectory, q = ', num2str(q)];
    title(Str, 'FontSize', 12, 'FontWeight', 'bold');    
    xlabel('X position in Cartesian Coordinate', 'FontSize', 10, ...
        'FontWeight','bold');
    ylabel('Y position in Cartesian Coordinate', 'FontSize', 10, ...
        'FontWeight','bold'); 
    
figure,
plot(T, ones(size(T))*10); hold on;
plot(T, (X_EKF(3,:).^2 + X_EKF(4,:).^2).^0.5, ':', 'color', 'red', ...
    'LineWidth', 2 );
    legend('Reference Velocity', 'EKF model velocuty', 'Location', 'Best');
    title('The EKF model velocity', 'FontSize', 12, 'FontWeight', 'bold');    
    xlabel('Time (sec)', 'FontSize', 10, 'FontWeight','bold');
    ylabel('Velocity (m/s)', 'FontSize', 10, 'FontWeight','bold');  
