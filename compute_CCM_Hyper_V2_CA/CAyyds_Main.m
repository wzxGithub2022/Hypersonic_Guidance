%基于在线算法框架，对CA进行重构
%寄，非线性的点撒下去，msspoly根本找不到收缩的CCM，还是得限定到一个小区域内，比线性化也差别不大
clear;clc;
warning off
close all

% ref是参考轨迹，x是射程，y是高度
global t_r u_r y_r h_r v_r gamma_r feedback disturb
global R0 g
R0 = 10*10^3;
g = 9.81; 

load('Trajectory_normalization.mat')
t_r = Trajectory_normalization(:,1);
y_r = Trajectory_normalization(:,2);
h_r = Trajectory_normalization(:,3);
v_r = Trajectory_normalization(:,4);
gamma_r = Trajectory_normalization(:,5);
u_r = Trajectory_normalization(:,6);

feedback = 1;   %是否加控制的总开关
disturb = 1;    %是否加干扰的总开关

% 控制Maxstep和RelTol来保证较少的计算点，提高计算速度
options = odeset('Maxstep', 0.1, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %返回了这些点上的值
[T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);

y_fig = X(:,1)*R0/1000;     %km
h_fig = X(:,2)*R0/1000;     %km
v_fig = X(:,3)*sqrt(R0*g);
gama_fig = rad2deg(X(:,4));

y_std_fig = y_r*R0/1000;     %km
h_std_fig = h_r*R0/1000;     %km
v_std_fig = v_r*sqrt(R0*g);     %m/s
gama_std_fig = rad2deg(gamma_r);    %度

figure(1)
plot(y_fig,h_fig,'LineWidth',1)
hold on
plot(y_std_fig,h_std_fig,'--','Color',[0.89,0.09,0.05],'LineWidth',1.5)
xlabel('relative lateral distance/(km)')
ylabel('height/(km)')
ylim([0 20.5])
legend('CCM tracking','nominal')
title('二维轨迹跟踪')

disp('射程偏差km');
disp(y_fig(end));
disp('高度偏差km');
disp(h_fig(end));
disp('落速m/s');
disp(v_fig(end));
disp('落角/度');
disp(gama_fig(end));