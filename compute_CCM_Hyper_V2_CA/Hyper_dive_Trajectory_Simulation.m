%%
%-------------------------------------------------------------------------%
% normalization里有所有归一化信息，t,y,h,v,gama
% M阵是拿归一化动力学算出来的，但是你仿真理应拿真是动力学来算
% 没有强力支持你归一化算出的CCM是能直接使用在原始动力学的，但你确实可以通过加密再解密的过程得到实际控制量
%-------------------------------------------------------------------------%
clc;clear;close all;
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor
global t0 t0_disp %为了观察程序运行到什么地步

R0 = 10*10^3;
g = 9.81; 

load('Trajectory_normalization.mat')
load('CCM_upper_CA.mat')
t_nor = Trajectory_normalization(:,1);                %控制时间序列
y_nor = Trajectory_normalization(:,2);
h_nor = Trajectory_normalization(:,3);
v_nor = Trajectory_normalization(:,4);
gama_nor = Trajectory_normalization(:,5);
alpha_nor = Trajectory_normalization(:,6);                %控制量u
CCM_nor = CCM_upper;

% 归一化的初始状态
y1_initial = Trajectory_normalization(1,2);
h1_initial = Trajectory_normalization(1,3);
v1_initial = Trajectory_normalization(1,4);
gama1_initial = Trajectory_normalization(1,5);
alpha1_initial = Trajectory_normalization(1,6);
state1_initial = [y1_initial, h1_initial, v1_initial, gama1_initial, alpha1_initial];         %初值航程35km 高度20km 速度 航迹角（弹道倾角)

tspan = [0 t_nor(end)];                            %时间区间
global step
step = 0.005/(sqrt(R0/g));                       %步长

options = odeset('RelTol',1e-9,'MaxStep',step);
[t,state1] = ode45(@Hyper_dive_Dynamitics_2D_function,tspan,state1_initial,options);  %二维动力学写在这里了
%% 跑一跑，验证CA的RelTol是否合适，控制是否真的大出天际
% CA能跑RelTol的1e-9，且相比于(157.8m,167.9m)误差有减小，控制也没有大出天际
% 射程偏差km
%     0.1111
% 高度偏差km
%    -0.1416
%%
%为了绘图，反归一化
t_fig = t*sqrt(R0/g);       %单位：s
y_fig = state1(:,1)/1000*R0;     %单位：km
h_fig = state1(:,2)/1000*R0;     %单位：km
v_fig = state1(:,3)*sqrt(R0*g);  %单位：m/s
gama_fig = rad2deg(state1(:,4)); %单位：度
alpha_fig = rad2deg(state1(:,5)); %单位：度

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

subplot(222),
plot(t_fig,alpha_fig,'LineWidth',2);
xlabel('时间t/s');ylabel('攻角\alpha/度');
title('攻角随时间的变化规律');
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
disp('射程偏差km');
disp(y_fig(end));
disp('高度偏差km');
disp(h_fig(end));