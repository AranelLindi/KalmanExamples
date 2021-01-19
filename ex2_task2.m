clear all
close all
load lrl_ex

%% Load sample data
% y_lidar       raw lidar measurements [cm]
% time_lidar    timestamps for lidar measurements
% y_accel       raw accelerometer measurements [m/s]
% time_accel    timestamps for accelerometer measurements
% 
% the following values are just for reference! 
% true_z        true height
% true_vz       true vertical speed
% true_accel    true vertical acceleration

time_accel = 0:0.005:(length(y_accel)-1)*0.005;
time_lidar = 0:0.05:(length(y_lidar)-1)*0.05;

% Calculate accel sample time
ST_accel = mean(diff(time_accel));

% Zwischenarray um Lidar-Graphen auf 100 Sekunden zu strecken:
lidar_func = zeros(1, length(y_accel));
last_lidar = 0;

%% Initial state
x = zeros(2,length(y_accel));
x(:,1) = [0; 0];

%% Initial system covariance 
P = [10 0; 0 10];

%% The system state propagation matrix F
dt = ST_accel;
F = [1 dt; 0 1];

%% The system input matrix G
G = [1/2*dt*dt; dt];

%% The measurement prediction matrix H
H = [100 0];

%% Measurement noise lidar
R_lidar = var(y_lidar(1:200));

%% Measurement noise accel
R_accel = var(y_accel(1:2000));

%% System prediction noise
Q = G * R_accel * G';

% iterate over time and measurments 
% il = index lidar, ia = index accel --> discrete time step
il = 1;
for ia = 1:length(y_accel)-1
    tic;
    %% Propagation accodring to system model
    % System state prediciton
    x(:,ia+1) = F * x(:, ia) + G * (y_accel(ia)-9.81); % u(k) = y_accel(k)
    % System state covariance prediction
    P = F * P * transpose(F) + Q;
    
    %% LiDar Measurement update
    if(time_accel(ia) >= time_lidar(il))
        
       % Measurement Prediciton 
       z = H * x(:, ia+1);
       % Measurement Residual
       v = y_lidar(il) - z;
       
       % Innovation covariance
       S = (H * P * transpose(H)) + R_lidar;
       % Filter Gain
       W = P * transpose(H) * inv(S);

       % Updated state estimate
       x(:,ia+1) = x(:, ia+1) + (W * v);
       % Updated state covariance
       P = (eye(2) - (W * H)) * P;
       
       % Increment lidar index
       il = il + 1;
       
       last_lidar = y_lidar(il)/100;
    end
    lidar_func(ia) = last_lidar;
    rt(ia) = toc;
end
lidar_func(end) = last_lidar;

disp(['Mean Runtime = ' num2str(mean(rt))])
disp(['MSE h = ' num2str(mean((x(1,:)'-true_z).^2))])
disp(['MSE v = ' num2str(mean((x(2,:)'-true_vz).^2))])

%% Visualize Results
% Plot the position
figure(1);
plot(x(1,:), 'LineWidth', 1); % Plot: geschätzter Zustand
hold on
plot(true_z, 'LineWidth', 2); % Plot: wahrer Zustand
val_lidar = y_lidar/100;
plot(lidar_func, 'LineWidth', 1); % Plot: Lidar-Werte
hold off
title('Verlauf der Höhe');
legend('geschätzter Systemzustand', 'wahrer Wert', 'Messwerte Abstandslaser');
legend('location', 'southeast');
% Hier: X-Achse anpassen, sodass sie anstatt Anzahl der Messungen, Zeitraum
% abbildet. Dazu mit Sampletime Beschleunigungssensor (ST_accel)
% multiplizieren:
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt*ST_accel);

xlabel('Zeit [s]');
ylabel('Höhe [m]');

%saveas(figure(1), 'Position.eps');
%clf;

% ##################################################

% Plot the speed
figure(2);
plot(x(2,:), 'LineWidth', 1); % Plot: geschätzter Systemzustand
hold on
plot(true_vz, 'LineWidth', 2); % Plot: wahrer Wert
hold off
title('Verlauf der Geschwindigkeit');
xlabel('Zeit [s]');
ylabel('Geschwindigkeit [m/s]');
legend('geschätzter Systemzustand', 'wahrer Wert');
% Hier: X-Achse anpassen, sodass sie anstatt Anzahl der Messungen, Zeitraum
% abbildet. Dazu mit Sampletime Beschleunigungssensor (ST_accel)
% multiplizieren:
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt*ST_accel);

%saveas(figure(2), 'Speed.eps');
%clf;

% ##################################################

% Plot the acceleration
figure(3);
y_accel_offset = y_accel(:,1) - 9.81; % Hier "Offset" (g) abziehen um schönere Ausgabe zu erhalten
plot(y_accel_offset, 'LineWidth', 1); % Plot: Messwerte Beschleunigungsssensor
hold on
plot(x(2,:), 'LineWidth', 1); % Plot: geschätzter Systemzustand
plot(true_accel, 'LineWidth', 2); % Plot: wahrer Wert
hold off
title('Verlauf der Beschleunigung');
xlabel('Zeit [s]');
ylabel('Beschleunigung [m/s^2]');
legend('Messwerte Beschleunigungssensor (abzgl. g)', 'geschätzter Systemzustand', 'wahrer Wert');
legend('location', 'southeast');
% Hier: X-Achse anpassen, sodass sie anstatt Anzahl der Messungen, Zeitraum
% abbildet. Dazu mit Sampletime Beschleunigungssensor (ST_accel)
% multiplizieren:
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt*ST_accel);

%saveas(figure(3), 'Accel.eps');
%clf;

