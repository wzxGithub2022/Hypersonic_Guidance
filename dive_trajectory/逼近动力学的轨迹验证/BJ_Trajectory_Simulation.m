%%
%-------------------------------------------------------------------------%
% normalization里有所有归一化信息，t,y,h,v,gama
% M阵是拿归一化动力学算出来的，但是你仿真理应拿真是动力学来算
% 没有强力支持你归一化算出的CCM是能直接使用在原始动力学的，但你确实可以通过加密再解密的过程得到实际控制量
%-------------------------------------------------------------------------%
clc;clear;close all;
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor
global t0 t0_disp %为了观察程序运行到什么地步

R0 = 10*10^3;
g = 9.81; 

load('Trajectory_normalization.mat')
t_nor = Trajectory_normalization(:,1);                %控制时间序列
y_nor = Trajectory_normalization(:,2);
h_nor = Trajectory_normalization(:,3);
v_nor = Trajectory_normalization(:,4);
gama_nor = Trajectory_normalization(:,5);
alpha_nor = Trajectory_normalization(:,6);                %控制量u

% 归一化的初始状态
y1_initial = Trajectory_normalization(1,2);
h1_initial = Trajectory_normalization(1,3);
v1_initial = Trajectory_normalization(1,4);
gama1_initial = Trajectory_normalization(1,5);
alpha1_initial = Trajectory_normalization(1,6);
state1_initial = [y1_initial, h1_initial, v1_initial, gama1_initial];         %初值航程35km 高度20km 速度 航迹角（弹道倾角)

%时间区间，看起来dt和dt1的问题还需要甄别
%你就是拿t_nor在仿真的，为什么要在程序中部分地恢复真实时间，已经换了一个世界了
% tspan = [ 0,  t_nor(end) * (sqrt(R0/g) ];
tspan = [ 0,  t_nor(end) ];

global step
step = 0.005;                       %步长

options = odeset('RelTol',1e-6,'MaxStep',step);
[t,state1] = ode45(@BJ_Dynamitics_2D_function,tspan,state1_initial,options);  %二维动力学写在这里了

%%
%为了绘图，反归一化
t_fig = t*sqrt(R0/g);       %单位：s
y_fig = state1(:,1)/1000*R0;     %单位：km
h_fig = state1(:,2)/1000*R0;     %单位：km
v_fig = state1(:,3)*sqrt(R0*g);  %单位：m/s
gama_fig = rad2deg(state1(:,4)); %单位：度

figure(1),
plot(t_fig,h_fig,'-ob',t_fig,y_fig,'-*r');
xlabel('时间t/s');ylabel('高度与射程/km');legend('高度h','射程y');
title('高度与射程随时间的变化曲线');
grid on

figure(2),subplot(221),
plot(y_fig,h_fig,'linewidth',2);
xlabel('射程y/km');ylabel('高度h/km');
title('二维平面最优轨迹');
grid on

subplot(223),
plot(t_fig,v_fig,'linewidth',3);
xlabel('时间t/s');ylabel('速度v/(m/s)');
title('速度随时间的变化曲线');
grid on

subplot(224),
plot(t_fig,gama_fig,'linewidth',3);
xlabel('时间t/s');ylabel('弹道倾角/度');
title('弹道倾角随时间的变化曲线');
grid on

%%
%为了绘图，反归一化
t_norf = t_nor*sqrt(R0/g);       %单位：s
y_norf = y_nor/1000*R0;     %单位：km
h_norf = h_nor/1000*R0;     %单位：km
v_norf = v_nor*sqrt(R0*g);  %单位：m/s
gama_norf = rad2deg(gama_nor); %单位：度

figure(3),
subplot(221),
plot(t_norf,y_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,y_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');ylabel('射程y/km');
legend('标称轨迹','逼近动力学仿真轨迹');
title('逼近动力学效果图');
grid on;

subplot(222),
plot(t_norf,h_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');ylabel('高度h/km');
legend('标称轨迹','逼近动力学仿真轨迹');
title('逼近动力学效果图');
grid on;

subplot(223),
plot(t_norf,v_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,v_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');ylabel('速度v/(m/s)');
legend('标称轨迹','逼近动力学仿真轨迹');
title('逼近动力学效果图');
grid on;

subplot(224),
plot(t_norf,gama_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,gama_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');ylabel('弹道倾角/度');
legend('标称轨迹','逼近动力学仿真轨迹');
title('逼近动力学效果图');
grid on;

%%
disp('射程偏差km');
disp(y_fig(end));
disp('高度偏差km');
disp(h_fig(end));
disp('落速m/s');
disp(v_fig(end));
disp('落角/度');
disp(gama_fig(end));
%%
disp('射程偏差km');
disp(y_norf(end));
disp('高度偏差km');
disp(h_norf(end));
disp('落速m/s');
disp(v_norf(end));
disp('落角/度');
disp(gama_norf(end));
%% 用逼近动力学，这成了印度布朗弹了
% 射程偏差km
%     4.4208
% 高度偏差km
%    -1.4095
%% 狗屁，这精度还可以接受，虽然单点算dot误差最大能到4%，但整体积分拉回来了一些
% 看看这个落点，现在对制导和导引的差别有点了啊
% 射程偏差km
%     0.1386
% 高度偏差km
%    -0.0220