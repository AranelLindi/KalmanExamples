%% Load data
clear all;
load data;
% Calculate difference between time steps
dt = diff(time);
mean_dt = mean(dt);

% figure(1);
% subplot(3,1,1);
% plot(time, acc);
% grid on;
% xlabel('Zeit [s]');
% ylabel('Beschleunigung [g]');
% title('Messwerte Beschleunigung');
% legend('x-Achse', 'y-Achse', 'z-Achse');
% legend('location', 'northwest');
% 
% subplot(3,1,2);
% plot(time, gyro);
% grid on;
% xlabel('Zeit [s]');
% ylabel('Drehrate [rad/s]');
% title('Messwerte Gyroskop');
% legend('x-Achse', 'y-Achse', 'z-Achse');
% legend('location', 'southwest');
% 
% subplot(3,1,3);
% plot(time, mag);
% grid on;
% xlabel('Zeit [s]');
% ylabel('mag. Flussdichte [gauss]');
% title('Messwerte Magnetometer');
% legend('x-Achse', 'y-Achse', 'z-Achse');
% legend('location', 'southeast');

%% Calculate initial condition
x = zeros(10,length(time));

% Samples used to calculate initial condition
% Assume platform is not moving for num_init samples
num_init = 150;	% 1s is enough
% 150 Datensätze entsprechen ungefähr 1 Sek.

% Mean values for num_init measurements
mean_gyro = mean(gyro(:,1:num_init),2);     % 3x1
mean_acc = -mean(acc(:,1:num_init),2);       % 3x1 % Achtung! Minus hinzugefügt! (siehe Folien S. 36)
mean_mag = mean(mag(:,1:num_init),2);       % 3x1

% Calculate initial roll angle
roll_init = atan2(mean_acc(2), mean_acc(3));

% Calculate initial pitch angle
pitch_init = atan2(-mean_acc(1), sqrt(mean_acc(2)^2+mean_acc(3)^2));

% Calculate initial yaw angle
cp = cos(pitch_init);
sp = sin(pitch_init);
cr = cos(roll_init);
sr = sin(roll_init);
% Rotation matrix: roll around x, pitch around y
R = [cp sp*sr cp*sr; 0 cr -sr; -sp cp*sr cp*cr];
mw = R * mean_mag;
yaw_init = atan2(-mw(2), mw(1));

% Calculate initial quaternion from yaw pitch and roll
q_init = euler2quat(yaw_init, pitch_init, roll_init);

% Set initial state
x(:,1) = [q_init(1); q_init(2); q_init(3); q_init(4); 0; 0; 0; mean_gyro(1); mean_gyro(2); mean_gyro(3)];         % 10x1

% Initial System State Covariance
P = 10*eye(10);                               % 10x10

%% Noise caracteristics
% System noise
v_bias = 1.0e-11;       % variance of bias
lambda = 1.0e-3;        % time coefficient for bias change

% Sensor noise

% (Es wird ein Zeitraum gewählt, der ungefähr im ersten Drittel während des
% Fluges stattfindet)
s_gyro = std(gyro(7625:10000));        % std of gyro

s_acc = std(acc(7625:10000));
mean_yaw = atan2(-mag(2,7625:10000), mag(1,7625:10000));

s_yaw  = std(mean_yaw);        % std of Magnetometer yaw

%% Measurement Noise
% Noise matrix of Magnetometer Yaw mesurement
Ryaw = s_yaw^2;                             % 1x1
% Noise matrix of Magnetometer Acc mesurement
Racc = diag([s_acc^2, s_acc^2, s_acc^2]);   % 3x3
% Noise matrix R
R =  [Racc , zeros(3,1); ...                   % 4x4
      zeros(1,3), Ryaw];

%% System Noise Q
% Projection matrix Fu
Fu = [0 0 0 0 0 0;...
      0 0 0 0 0 0;...
      0 0 0 0 0 0;...
      0 0 0 0 0 0;...
      1 0 0 0 0 0; ...
      0 1 0 0 0 0; ...
      0 0 1 0 0 0; ...
      0 0 0 1 0 0; ...
      0 0 0 0 1 0; ...
      0 0 0 0 0 1];              %10x6

% System input noise
U = diag([s_gyro^2, s_gyro^2, s_gyro^2, v_bias, v_bias, v_bias]); % 6x6
% Noise matrix Q
Q = Fu * U * Fu';               

%% Quaternion EKF
for k=1:length(time)-1
    %% State prediction
    % Auxiliary variables
    % Werte aus der k-Iteration werden benötigt (siehe unten) um neue Werte
    % für k+1 zu berechnen:
    q0 = x(1,k);
    q1 = x(2,k);
    q2 = x(3,k);
    q3 = x(4,k);
    wx = x(5,k);
    wy = x(6,k);
    wz = x(7,k);
    xgx = x(8,k);
    xgy = x(9,k);
    xgz = x(10,k);
    
    % Non-linear system state prediciton x_(k+1) = f(x_(k), gyro)
    x(1,k+1)  = q0 + 1/2 * (-(q1*wx+q2*wy+q3*wz))*dt(k);
    x(2,k+1)  = q1 + 1/2 * (q0*wx-q3*wy+q2*wz)*dt(k);
    x(3,k+1)  = q2 + 1/2 * (q3*wx+q0*wy-q1*wz)*dt(k);
    x(4,k+1)  = q3 + 1/2 * (-q2*wx+q1*wy+q0*wz)*dt(k);
    x(5,k+1)  = gyro(1,k)-xgx; % Neue Messwerte von Gyro abholen und Bias abziehen
    x(6,k+1)  = gyro(2,k)-xgy; % Wichtig! Da hier Messwerte abgeholt werden, muss in Fx der Systemeingang nicht gesetzt werden
    x(7,k+1)  = gyro(3,k)-xgz; % "   " 
    x(8,k+1)  = xgx; % Bleibt konstant
    x(9,k+1)  = xgy; % ""
    x(10,k+1) = xgz; % ""
    
    % Normalize quaternion to compensate numerical inaccuracies
    x(1:4,k+1) = x(1:4,k+1) ./ sqrt(sum(x(1:4,k+1).^2));
    
    % Jacobian of f(x,u): Fx   % 10x10
    Fx = [1 0 0 0 (-1/2*q1*dt(k)) (-1/2*q2*dt(k)) (-1/2*q3*dt(k)) 0 0 0; ...
          0 1 0 0 (1/2*q0*dt(k)) (-1/2*q3*dt(k)) (1/2*q2*dt(k)) 0 0 0; ...
          0 0 1 0 (1/2*q3*dt(k)) (1/2*q0*dt(k)) (-1/2*q1*dt(k)) 0 0 0; ...
          0 0 0 1 (-1/2*q2*dt(k)) (1/2*q1*dt(k)) (1/2*q0*dt(k)) 0 0 0; ...
          0 0 0 0 0 0 0 -1 0 0; ...
          0 0 0 0 0 0 0 0 -1 0; ...
          0 0 0 0 0 0 0 0 0 -1; ...
          0 0 0 0 0 0 0 1 0 0; ...
          0 0 0 0 0 0 0 0 1 0; ...
          0 0 0 0 0 0 0 0 0 1];
     
    % System state covariance prediction
    P = Fx * P * Fx' + Q; 
    
    % Update auxiliary variables
    q0 = x(1,k+1);
    q1 = x(2,k+1);
    q2 = x(3,k+1);
    q3 = x(4,k+1);
    
    %% Measurement Update
    % Rotationmatrix from body frame to navigation frame
    R_bn = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3), 2*(q0*q2+q1*q3); ...
            2*(q1*q2+q0*q3), (q0^2-q1^2+q2^2-q3^2), 2*(q2*q3-q0*q1); ...
            2*(q1*q3-q0*q2), 2*(q0*q1+q2*q3), (q0^2-q1^2-q2^2+q3^2)]; 
        
    % Rotationmatrix from navigation to body frame
    R_nb = R_bn';

    %% Accelerometer
    % Accelerometer measurement yacc
    yacc = acc(:,k);
    % Non-linear measurement prediciton of accelerometer zacc
    % Estimate gravitational vector as measured in the body frame
    zacc = R_nb * yacc;

    % Jacobian of zacc = hacc(x): Hacc      %3x10
    Hacc = [2*q2 -2*q3 2*q0 -2*q1 0 0 0 0 0 0; ...
            -2*q1 -2*q0 -2*q3 -2*q2 0 0 0 0 0 0; ...
            -2*q0 2*q1 2*q2 -2*q3 0 0 0 0 0 0];

    %% Magnetometer
    % Magnetometer measurement ymag
    ymag = mag(:,k);
    % Rotate Magnetometer measurement into navigation frame
    mn = R_bn * ymag;
    % Project mn to x-y-plane
    mnh = [mn(1); mn(2); 0];
    % Rotate projection back to body frame
    mb = R_nb * mnh;
    % Calculate yaw measurement yyaw from new vector mb 
    yyaw = atan2(-mb(2), mb(1));
    
    % Non-linear measurement prediciton of yaw angle zyaw
    % Estimate gravitational vector as measured in the body frame
    zyaw = atan2(2*(q0*q3+q1*q2), (1-2*(q2^2+q3^2)));  
    
    % Jacobian of zyaw = hyaw(x): Hyaw          %1x10
    Hyaw = zeros(1,10);
    Hyaw(1,1) = (2*q3*(1-2*(q2^2+q3^2)))/(4*(q0*q3+q1*q2)^2+(1-2*(q2^2+q3^2))^2);
    Hyaw(1,2) = (2*q2*(1-2*(q2^2+q3^2)))/(4*(q0*q3+q1*q2)^2+(1-2*(q2^2+q3^2))^2);
    Hyaw(1,3) = (8*q2*(q0*q3+q1*q2)+2*q1*(1-2*(q2^2+q3^2)))/(4*(q0*q3+q1*q2)^2+(1-2*(q2^2+q3^2))^2);
    Hyaw(1,4) = (8*q3*(q0*q3+q1*q2)+2*q0*(1-2*(q2^2+q3^2)))/(4*(q0*q3+q1*q2)^2+(1-2*(q2^2+q3^2))^2);

    %% Common update
    %% Combined update
    % Combined Measurement Prediction
    z = [zacc; zyaw];
    % Combined Measurement Vector
    y = [yacc; yyaw];
    % Combined Measurement Jacobian
    H = [Hacc; Hyaw];

    % Combined Measurement Residual
    v = (y-z);

    % Limit yaw error to be between -180° and 180°
    if (v(4) > pi)
    v(4) = v(4) - 2*pi;
    elseif (v(4) < -pi)
     v(4) = v(4) + 2*pi;
    end
    
    % Innovation covariance
    S = H * P * H' + R;
    
    % Filter Gain
    W = P * H' * inv(S);
    
    % Updated state estimate
    x(:,k+1) = x(:,k+1) + (W * v);
    % Normalize quaternion to compensate numerical inaccuracies
    x(1:4,k+1) = x(1:4,k+1) ./ sqrt(sum(x(1:4,k+1).^2));
    
    % Updated state covariance
    P = (eye(10)-W*H)*P;
end

%% Plot quaternion as EULER angles
% Get the quaternion from state
quat = x(1:4,:)';
% Convert to euler angles
[y_ekf, p_ekf, r_ekf] = quat2euler(quat);
[y_true, p_true, r_true] = quat2euler(true_q);

figure(2);
subplot(3,1,1);
plot(time, [y_ekf, p_ekf, r_ekf], 'LineWidth', 2)
grid on;
xlabel('Zeit [s]');
ylabel('YPR [°]');
title('Geschätzter Systemzustand durch QEKF');
legend('Yaw', 'Pitch', 'Roll');

subplot(3,1,2);
plot(time, [y_true, p_true, r_true], 'LineWidth', 2)
grid on;
xlabel('Zeit [s]');
ylabel('YPR [°]');
title('Wahrer Wert');
legend('Yaw', 'Pitch', 'Roll');

subplot(3,1,3);
plot(time, [y_ekf, p_ekf, r_ekf] - [y_true, p_true, r_true], 'LineWidth', 2)
grid on;
xlabel('Zeit [s]');
ylabel('YPR [°]');
title('Fehler/Abweichung');
legend('Yaw', 'Pitch', 'Roll');
legend('location', 'southeast');

