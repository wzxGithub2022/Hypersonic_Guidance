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
%%
feedback = 1;   %是否加控制的总开关
disturb = 1;    %是否加干扰的总开关

% 控制Maxstep和RelTol来保证较少的计算点，提高计算速度
% step = 0.005/sqrt(R0/g);
step = 0.1;
options = odeset('Maxstep', step, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %返回了这些点上的值
[T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);

t_fig = tspan*sqrt(R0/g);
y_fig = X(:,1)*R0/1000;     %km
h_fig = X(:,2)*R0/1000;     %km
v_fig = X(:,3)*sqrt(R0*g);
gama_fig = rad2deg(X(:,4));
%%
feedback = 0;   %是否加控制的总开关
disturb = 0;    %是否加干扰的总开关

% 控制Maxstep和RelTol来保证较少的计算点，提高计算速度
% step = 0.005/sqrt(R0/g);
step = 0.1;
options = odeset('Maxstep', step, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %返回了这些点上的值
[T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);

t_std_fig = tspan*sqrt(R0/g);
y_std_fig = X(:,1)*R0/1000;     %km
h_std_fig = X(:,2)*R0/1000;     %km
v_std_fig = X(:,3)*sqrt(R0*g);
gama_std_fig = rad2deg(X(:,4));
%%
feedback = 0;   %是否加控制的总开关
disturb = 1;    %是否加干扰的总开关

% 控制Maxstep和RelTol来保证较少的计算点，提高计算速度
% step = 0.005/sqrt(R0/g);
step = 0.1;
options = odeset('Maxstep', step, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %返回了这些点上的值
[T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);

t_dtb_fig = tspan*sqrt(R0/g);
y_dtb_fig = X(:,1)*R0/1000;     %km
h_dtb_fig = X(:,2)*R0/1000;     %km
v_dtb_fig = X(:,3)*sqrt(R0*g);
gama_dtb_fig = rad2deg(X(:,4));
%%
figure(1)
plot(y_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(y_dtb_fig,h_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(y_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('射程y/km');
ylabel('高度h/km');
axis([-1 35 0 21]);
legend('标称轨迹','干扰轨迹','制导轨迹');
title('轨迹仿真效果');
grid on

figure(2)
subplot(221)
plot(t_std_fig,y_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,y_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,y_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('射程y/km');
legend('标称轨迹','干扰轨迹','制导轨迹');
title('轨迹仿真效果');
grid on

subplot(222)
plot(t_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,h_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('高度h/km');
legend('标称轨迹','干扰轨迹','制导轨迹');
title('轨迹仿真效果');
grid on

subplot(223)
plot(t_std_fig,v_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,v_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,v_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('速度v/(m/s)');
legend('标称轨迹','干扰轨迹','制导轨迹');
title('轨迹仿真效果');
grid on

subplot(224)
plot(t_std_fig,gama_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,gama_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,gama_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('弹道倾角\gamma/度');
legend('标称轨迹','干扰轨迹','制导轨迹');
title('轨迹仿真效果');
grid on

%%
Bias = {'射程偏差km';'高度偏差km';'落速m/s';'落角/度'};
Norminal = [y_std_fig(end);h_std_fig(end);v_std_fig(end);gama_std_fig(end)];
Guidance = [y_fig(end);h_fig(end);v_fig(end);gama_fig(end)];
Disturb = [y_dtb_fig(end);h_dtb_fig(end);v_dtb_fig(end);gama_dtb_fig(end)];
T = table(Bias,Norminal,Guidance,Disturb)
%%
% disp('制导-射程偏差km');
% disp(y_fig(end));
% disp('制导-高度偏差km');
% disp(h_fig(end));
% disp('制导-落速m/s');
% disp(v_fig(end));
% disp('制导-落角/度');
% disp(gama_fig(end));
% 
% disp('扰动-射程偏差km');
% disp(y_dtb_fig(end));
% disp('扰动-高度偏差km');
% disp(h_dtb_fig(end));
% disp('扰动-落速m/s');
% disp(v_dtb_fig(end));
% disp('扰动-落角/度');
% disp(gama_dtb_fig(end));
%% old
% %%
% % figure(1)
% % plot(y_fig,h_fig,'LineWidth',1)
% % hold on
% % plot(y_std_fig,h_std_fig,'--','Color',[0.89,0.09,0.05],'LineWidth',1.5)
% % xlabel('relative lateral distance/(km)')
% % ylabel('height/(km)')
% % ylim([0 20.5])
% % legend('CCM tracking','nominal')
% % title('二维轨迹跟踪')
% %%
% % 射程偏差km
% %    -0.0037
% % 高度偏差km
% %    -0.0475
% % 落速m/s
% %    1.0801e+03
% % 落角/度
% %   -65.4138