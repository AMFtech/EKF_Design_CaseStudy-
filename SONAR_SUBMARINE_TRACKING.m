% KALMAN FILTER EXMPLE FOR TRACKING SINGLE CART 
% CONSTANT VELOCITY ASSUMPTION IS MADE
% Assumptions: Theta is in Radian, Rho is in meter, Time is in sec and the
% velocuty is assumed to be constant through the experiement....
 
clc
clear all;
% close all;

Ts = 0.1;           % sampling time, dt = 0.1 sec
load object_trajectory.mat
RealPos = x_traj_pos;
[m, n] = size(RealPos);

% Calculate the Real Velocity Components
PosDelayed = horzcat(RealPos(:,2:end), [0 ; 0]);
RealVel = ( PosDelayed - RealPos )/Ts;
% figure; 
% quiver(RealPos(1,:), RealPos(2,:), RealVel(1,:), RealVel(2,:));

% Process Covariance factor
global q;           
q = 0.10;

% Observation Covariance Matrix
global Sig_rho;
global Sig_theta;

Sig_rho   = 0.1;            % in meter
Sig_theta = 3 * (pi/180);   % in radian

ShowPlot = 0;

%% Generate 100 Data set (noisy trajectories)
% m-n-100 size matrix, each m-n matrix contains single data set
if (ShowPlot)
    figure;
end
[SonarDataSet, SonarDataSetCart] = SonarDataSetGen(ShowPlot);

%% EKF operation
[mZ, nZ, pZ] = size(SonarDataSet);
Sonar_i = zeros(mZ, nZ);

Xk_EKF  = zeros(4,n,pZ);
XY_MSE  = zeros(2,n);
Vxy_MSE = zeros(2,n);

Xk_MEAN = zeros(4,n);

Xk_vs_Zk = zeros(1,n);
Pk_Norm  = zeros(1,n);

for i = 1:100
    Sonar_i(:,:) = SonarDataSet(:,:,i);
    [X_EKFi, K_Factor_i, Pk_Norm_i] = myEKF(Sonar_i);
    Xk_EKF(:,:,i) = X_EKFi(:,:);    
    XY_MSE = XY_MSE + ( X_EKFi(1:2,:) -  RealPos ).^2;
    Vxy_MSE = Vxy_MSE + ( X_EKFi(3:4,:) -  RealVel ).^2;
    
    Xk_MEAN = Xk_MEAN + X_EKFi;
    Xk_vs_Zk = K_Factor_i + Xk_vs_Zk;
    Pk_Norm  = Pk_Norm + Pk_Norm_i;
end
    XY_RMSE     = sqrt(XY_MSE / 100);
    Vxy_RMSE    = sqrt(Vxy_MSE / 100);
    Xk_MEAN     = Xk_MEAN / 100;
    Xk_vs_Zk    = Xk_vs_Zk / 100;
    Pk_Norm     = Pk_Norm / 100;

%% Evaluate Root-Mean-Square-Error of Xk for 100 data set 

T = 0:Ts:Ts*(n-1);  % time course through the experiment

% Plot the Norm of Pk, the process Covariance Matrix ----------------------  
figure;
plot(T, Pk_Norm, '--b', 'LineWidth', 2 );
    Str = ['|P_k| in time over the trajectory, q = ' num2str(q)];
    title(Str,'FontSize', 12, ...
        'FontWeight','bold');
    xlabel('Time (sec)', 'FontSize', 11, ...
        'FontWeight','bold');
    ylabel('Factor', 'FontSize', 11, ...
        'FontWeight','bold');  

% Plot the ratio of norm( Hk*Pp*Hk') / norm(R) ----------------------------  
figure;
plot(T, Xk_vs_Zk, '--b', 'LineWidth', 2 );
    Str = ['The ratio of |H_k * P_p * H_k^t / |R|, q = ' num2str(q)];
    title(Str,'FontSize', 12, ...
        'FontWeight','bold');
    xlabel('Time (sec)', 'FontSize', 11, ...
        'FontWeight','bold');
    ylabel('Factor', 'FontSize', 11, ...
        'FontWeight','bold');  

% Plot mean EKF trajectory and mean EKF velocity --------------------------
figure;
subplot(2,1,1);
plot(RealPos(1,:), RealPos(2,:), '--r', 'LineWidth', 2 ); hold on;
plot(Xk_MEAN(1,:), Xk_MEAN(2,:),'b', 'LineWidth', 3 ); 
    legend('Real Trajectory', 'Mean Trajectory', 'Location', 'Best');
    Str = ['EKF Mean trajectory over 100 Data Set, q = ' num2str(q)];
    title(Str,'FontSize', 12, ...
        'FontWeight','bold');
    xlabel('x_k', 'FontSize', 11, ...
        'FontWeight','bold');
    ylabel('y_k', 'FontSize', 11, ...
        'FontWeight','bold');  
subplot(2,1,2);
plot(T(1:end-1), sqrt( Xk_MEAN(3,1:end-1).^2 + Xk_MEAN(4,1:end-1).^2 ), ...
    'r', 'LineWidth', 3 );
hold on; 
plot( T(1:end-1), 10*ones(1,n-1), '--r') ;
    legend('EKF mean velocuty','Reference velocity', 'Location', 'Best');
    Str = ['[x_k, y_k] Mean Square Error over 100 Data Set, q = ' num2str(q)];
    title(Str,'FontSize', 12, ...
        'FontWeight','bold');
    xlabel('Time (sec)', 'FontSize', 11, ...
        'FontWeight','bold');
    ylabel('Root Mean Square Error', 'FontSize', 11, ...
        'FontWeight','bold');    
%--------------------------------------------------------------------------


% Plot Root Mean Square Error of EKF position and velocity ----------------
figure;
subplot(2,1,1);
plot(T, XY_RMSE(1,:),'b', 'LineWidth', 3 );
hold on;
plot(T, XY_RMSE(2,:),'r', 'LineWidth', 3 );
legend('RMSE for x_k','RMSE for y_k', 'Location', 'Best');
    Str = ['Position Root Mean Square Error over 100 Data Set, q = ' num2str(q)];
    title(Str,'FontSize', 12, ...
        'FontWeight','bold');
    xlabel('Time (sec)', 'FontSize', 11, ...
        'FontWeight','bold');
    ylabel('Root Mean Square Error', 'FontSize', 11, ...
        'FontWeight','bold');    
    
subplot(2,1,2);
plot(T(1:end-1), Vxy_RMSE(1,1:end-1),'b', 'LineWidth', 3 );
hold on;
plot(T(1:end-1), Vxy_RMSE(2,1:end-1),'r', 'LineWidth', 3 );
    legend('RMSE for V_x','RMSE for V_y', 'Location', 'Best');
    Str = ['Velocity Root Mean Square Error over 100 Data Set, q = ' num2str(q)];
    title(Str,'FontSize', 12, ...
        'FontWeight','bold');
    xlabel('Time (sec)', 'FontSize', 11, ...
        'FontWeight','bold');
    ylabel('Root Mean Square Error', 'FontSize', 11, ...
        'FontWeight','bold');    
%--------------------------------------------------------------------------

