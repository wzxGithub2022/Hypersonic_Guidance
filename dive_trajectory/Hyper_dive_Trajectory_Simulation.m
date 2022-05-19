%-------------------------------------------------------------------------%
%------------------ plot the trajectory with ode -------------------------%
%------------------ 2D_plain(y,z)-----------------------------------------%
%------------------ target location  (0,0) -------------------------------%
%------------------ initial location (y0,z0) -----------------------------%
%-------------------------------------------------------------------------%
clc;clear;close all;
%%
global t_u u

R0 = 10*10^3;
g0 = 9.81; 

load('Trajectory_control_information.mat')
t_u = time_control_sequence(:,1);                %控制时间序列
u   = time_control_sequence(:,2);                %控制量u

impact_angle = deg2rad(-65);                     %落角约束
tspan = [0 t_u(end)];                            %时间区间
y0 = [35*10^3 20*10^3 1750 deg2rad(-5)];         %初值航程35km 高度20km 速度 航迹角（弹道倾角)
step = 0.005;                                    %步长
options = odeset('RelTol',1e-6,'MaxStep',step);
[t,y] = ode45(@Hyper_dive_Dynamitics_2D_function,tspan,y0,options);  %二维动力学写在这里了

%%
%作图
y_sim = y(:,1)/1000;
h_sim = y(:,2)/1000;

load('Trajectory_opt_information.mat')
y_opt = trajectory_location_sequence(:,1);
h_opt = trajectory_location_sequence(:,2);

figure(1)
hold on,plot(y_opt,h_opt,'-*','Color',[0.85 0.325 0.098],'LineWidth',1);
hold on,plot(y_sim,h_sim,'--','Color',[0 0.447 0.741],'LineWidth',2);
xlabel('航程y/km');ylabel('高度h/km');
legend('优化解算轨迹','动力学仿真轨迹');
title('动力学仿真与优化解算轨迹对比');
grid on;

%%
%误差统计
fprintf('终端航程：%d (m)\n',  y(end,1));
fprintf('终端高度：%d (m)\n',  y(end,2));
fprintf('终端落速：%d (m/s)\n',y(end,3));
fprintf('落速偏差：%d (m/s)\n',y(end,3)-1080);
fprintf('终端落角：%d (度)\n',rad2deg(y(end,4)));
fprintf('落角偏差：%d (度)\n',rad2deg(y(end,4))-(-65));
%%
% 终端航程：-7.684552e-02 (m)
% 终端高度：3.606829e-02 (m)
% 终端落速：1.076287e+03 (m/s)
% 终端落角：-6.499940e+01 (度)