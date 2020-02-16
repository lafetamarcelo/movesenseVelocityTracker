
%%
filename = ".\Out of water.csv";
fileid = fopen(filename);
c = textscan(fileid,'%s %s %s %s %s %s %s','Delimiter',',','EmptyValue',0);
fclose(fileid);
clc;

%%
X = data.Xms2(2:end) + "." + data.Yms2(2:end);
Y = data.Zms2(2:end) + "." + data.VarName5(2:end);
Z = data.VarName6(2:end) + "." + data.VarName7(2:end);
t = data.Timestampms(2:end);

%%
X = str2num(char(X));
Y = str2num(char(Y));
Z = str2num(char(Z));
t = str2num(char(t));
t = (t - t(1))/1000;

figure(2)
subplot(2,1,1)
plot(t,X);
hold on;
plot(t,Y);
plot(t,Z);

hold off;
subplot(2,1,2);
plot(diff(t));

