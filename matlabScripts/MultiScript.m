clear; clc;

%% Read the multiple sensor data from csv file 

filename = "..\database\RotationTest.csv";
fileid = fopen(filename);
data = textscan(fileid,'%s %s %s %s %s %s %s %s',...
    'Delimiter',',', ...
    'EmptyValue',0);
fclose(fileid);


%%

sensor_data_size = length(c{2}) - 1;

% Create the linear acceleration 
liAcc.t = 0; liAcc.x = 0;
liAcc.y = 0; liAcc.z = 0;
liAcc.size = 1;

anVel.t = 0; anVel.x = 0;
anVel.y = 0; anVel.z = 0;
anVel.size = 1;

for i = 2 : (sensor_data_size + 1)
    if strcmp(c{1}(i),'LinearAcc')
       liAcc.t(liAcc.size) = str2num(string(c{2}(i)));
       liAcc.x(liAcc.size) = str2num(string(c{3}(i)) + "." + string(c{4}(i)));
       liAcc.y(liAcc.size) = str2num(string(c{5}(i)) + "." + string(c{6}(i)));
       liAcc.z(liAcc.size) = str2num(string(c{7}(i)) + "." + string(c{8}(i)));
       liAcc.size = liAcc.size + 1;
    else
       anVel.t(anVel.size) = str2num(string(c{2}(i)));
       anVel.x(anVel.size) = str2num(string(c{3}(i)) + "." + string(c{4}(i)));
       anVel.y(anVel.size) = str2num(string(c{5}(i)) + "." + string(c{6}(i)));
       anVel.z(anVel.size) = str2num(string(c{7}(i)) + "." + string(c{8}(i)));
       anVel.size = anVel.size +1;
    end
end

%%

lw = 2.0;

figure(1)
subplot(2,1,1)
hold on;
plot(anVel.t,anVel.x, 'LineWidth', lw);
plot(anVel.t,anVel.y, 'LineWidth', lw);
plot(anVel.t,anVel.z, 'LineWidth', lw);
legend([{'X'},{'Y'},{'Z'}])

hold off;
subplot(2,1,2);
plot(diff(anVel.t), 'LineWidth', lw);

figure(2)
subplot(2,1,1)
hold on;
plot(liAcc.t,liAcc.x, 'LineWidth', lw);
plot(liAcc.t,liAcc.y, 'LineWidth', lw);
plot(liAcc.t,liAcc.z, 'LineWidth', lw);
legend([{'X'},{'Y'},{'Z'}])

hold off;
subplot(2,1,2);
plot(diff(liAcc.t), 'LineWidth', lw);

