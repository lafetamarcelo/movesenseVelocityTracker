

%% Read the sensors data
clear; close all; clc;

studyCaseFilename = '..\database\RotationTest.csv';
multipleSensors   = true; 

liAcc = sensor('LinearAcc', studyCaseFilename, multipleSensors);
anVel = sensor('AngularVelocity', studyCaseFilename, multipleSensors);

%% Manipulate the sensor information

liAcc.plot_serie('fignum', 1, ...
                 'ylabel', "Linear Acceleration (m/s^2)"); 
anVel.plot_serie('fignum', 11, ...
                 'ylabel', "Angular Velocity (degree/s)", ...
                 'legend', ["Pitch", "Roll", "Yaw"] );

liAcc.resample_series('Enhance', 5) % Resample raw data
anVel.resample_series('Enhance', 5) % Resample raw data

liAcc.plot_serie('fignum', 2, ...
                 'ylabel', "Linear Acceleration (m/s^2)"); 
anVel.plot_serie('fignum', 12, ...
                 'ylabel', "Angular Velocity (degree/s)", ...
                 'legend', ["Pitch", "Roll", "Yaw"] );

liAcc.denoise_series() % Denoise the resampled data
anVel.denoise_series() % Denoise the resampled data

liAcc.plot_serie('fignum', 3, ...
                 'ylabel', "Linear Acceleration (m/s^2)"); 
anVel.plot_serie('fignum', 13, ...
                 'ylabel', "Angular Velocity (degree/s)", ...
                 'legend', ["Pitch", "Roll", "Yaw"] );

             
%% Integrate the sensor signal and plot

signal = anVel.integrate();

figure(31); hold on;
plot(signal.t, signal.x, 'LineWidth',2);
plot(signal.t, signal.y, 'LineWidth',2);
plot(signal.t, signal.z, 'LineWidth',2);
legend(["Pitch","Plunge","Yaw"]);
xlabel("Time (s)"); ylabel("Degree");
set(gcf,'color','w');
grid on; hold off;



%% Plot and show the gravity vector components

% Create the acceleration data
acc_t = liAcc.Timeseries.t;
acc_b = [liAcc.Timeseries.x; ...
         liAcc.Timeseries.y; ...
         liAcc.Timeseries.z];

acc_i = zeros(3, length(acc_b));
for k = 1 : length(acc_b)
    
   roti_b = utils.rot_matrix(signal.x(k), signal.y(k), signal.z(k)); 
   acc_i(:,k) = roti_b * acc_b(:,k);
end

figure(32); hold on;
title("Inertial accelerations");
plot(acc_t, acc_i, 'LineWidth', 2);
plot(acc_t, 9.81*ones(1,length(acc_b)),'-.','LineWidth',2);
legend(["Inertial x","Inertial y","Inertial z", "Gravity Reference"]);
grid on; xlabel("Time (s)"); ylabel("m/s^2");
set(gcf,'color','w'); hold off;

%% Show the gravity vector moving on 3d

figure(33); hold on;
title("Gravity vector");
grid on; axis equal; view(145,20);
xl = [-5,5]; set(gcf,'color','w');
line(2*xl, [0,0], [0,0], 'LineWidth', 1.5, 'Color', 'k');
line([0,0], 2*xl, [0,0], 'LineWidth', 1.5, 'Color', 'k');
line([0,0], [0,0], 2*xl, 'LineWidth', 1.5, 'Color', 'k');

for k = 1 : length(acc_b)
     plot3(acc_i(1,k),acc_i(2,k),-acc_i(3,k),'.r');
     pause(0.05); 
end
hold off;



%% Plot and show the orientation moving

orient = zeros(3, length(acc_b));
for k = 1 : length(acc_b)    
   roti_b = utils.rot_matrix(signal.x(k), signal.y(k), signal.z(k)); 
   orient(:,k) = roti_b * [0; 0; 7];
end


figure(34); hold on;
title("Orientation vector");
grid on; axis equal; view(145,20);
xl = [-5,5]; set(gcf,'color','w');
line(2*xl, [0,0], [0,0], 'LineWidth', 1, 'Color', 'k');
line([0,0], 2*xl, [0,0], 'LineWidth', 1, 'Color', 'k');
line([0,0], [0,0], 2*xl, 'LineWidth', 1, 'Color', 'k');

for k = 1 : length(acc_b)
     plot3(orient(1,k),orient(2,k),-orient(3,k),'.r');
     pause(0.01); 
end
hold off;



%% 

serie_x = anVel.Timeseries.x;
serie_y = anVel.Timeseries.y;
serie_z = anVel.Timeseries.z;
serie_t = anVel.Timeseries.t;      
             
%%

%series_data = [anVel.Timeseries.x; ...
%               anVel.Timeseries.y; ...
%               anVel.Timeseries.z]';

series_data = [serie_x; ...
               serie_y; ...
               serie_z]';

ts_obj = iddata(series_data, [], anVel.sampleTime.Value/1000);
g = etfe(ts_obj);
res_data = zeros(3, length(g.Frequency));
size_array = [1, length(g.Frequency)];
res_data(1,:) = reshape(g.SpectrumData(1,1,:), size_array);
res_data(2,:) = reshape(g.SpectrumData(2,2,:), size_array);
res_data(3,:) = reshape(g.SpectrumData(3,3,:), size_array);

figure(100);
plot(g.Frequency, res_data, 'LineWidth',2);
%set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');


series_data = [anVel.Timeseries.x; ...
               anVel.Timeseries.y; ...
               anVel.Timeseries.z]';
ts_obj = iddata(series_data, [], anVel.sampleTime.Value/1000);
g2 = spa(ts_obj);
res_data = zeros(3, length(g2.Frequency));
size_array = [1, length(g2.Frequency)];
res_data(1,:) = reshape(g2.SpectrumData(1,1,:), size_array);
res_data(2,:) = reshape(g2.SpectrumData(2,2,:), size_array);
res_data(3,:) = reshape(g2.SpectrumData(3,3,:), size_array);

figure(200);
plot(g2.Frequency, res_data, 'LineWidth',2);
%set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

figure(300);
spectrum(g, g2);

figure(400);
h = bodeplot(g2);
showConfidence(h,3);

%% Compute the filtering of the signal


%% Compute the angles from the velocity


dt   = anVel.sampleTime.Value / 1000 ;
pos  = zeros(3,anVel.dataSize);
tpos = zeros(1,anVel.dataSize);
tpos(1) = anVel.Timeseries.t(1) / 1000;
for k = 1 : anVel.dataSize - 1
    
    pos(:,k+1) = pos(:,k) + dt * [anVel.Timeseries.x(k); ...
                                  anVel.Timeseries.y(k); ...
                                  anVel.Timeseries.z(k)];
    tpos(1,k+1) = k*dt + anVel.Timeseries.t(1)/1000;
end

figure(100);
title('Rotation Test -- P.P.Y.');
plot(tpos, pos, 'LineWidth', 2);
legend(["Pitch","Plunge","Yall"]);
ylabel('Degrees'); xlabel('Time (s)');

%% Compute the angles from the velocity


dt   = anVel.sampleTime.Value / 1000 ;
pos  = zeros(3,anVel.dataSize);
tpos = zeros(1,anVel.dataSize);
tpos(1) = anVel.Timeseries.t(1) / 1000;
for k = 1 : anVel.dataSize - 1
    
    pos(:,k+1) = pos(:,k) + dt * [serie_x(k); ...
                                  serie_y(k); ...
                                  serie_z(k)];
    tpos(1,k+1) = (k)*dt + anVel.Timeseries.t(1)/1000;
end

figure(100);
title('Rotation Test -- P.P.Y.');
plot(tpos, pos, 'LineWidth', 2);
legend(["Pitch","Plunge","Yall"]);
ylabel('Degrees'); xlabel('Time (s)');



%%
%%
%% One more resample

time = anVel.Timeseries.t;
data = [anVel.Timeseries.x; ...
        anVel.Timeseries.y; ...
        anVel.Timeseries.z]';
cur_freq = 1000 / (time(2) - time(1));

enhance = 5;
[rdata, rt] = resample(data,time, enhance * cur_freq);

%%

shape = length(rt);
phi   = zeros(3,shape); 
theta = zeros(3,shape); 
psi   = zeros(3,shape);
p_phi   = - pi/2;
p_theta = - pi;
p_psi   = 0;

dt = (rt(2) - rt(1)) / 1000;
for k = 1 : shape
    
    p_phi = [];
    
end


figure(100);
plot(rt, rdata, 'LineWidth', 2)

%% Compute the rotational matrix

% Sample trajectory
t = 0:0.1:10;
x = 5*cos(t);
y = 5*sin(t);
z = sin(t)+cos(t)+t;
plot3(x,y,z,'*r');
%if you want to animate it
for i=1:length(x)
  plot3(x(i),y(i),z(i),'*r');
  hold on;
  pause(0.01);
end