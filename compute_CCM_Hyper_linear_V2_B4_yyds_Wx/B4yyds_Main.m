%为了提高程序的可比性和运行效率，基于在线算法框架，对B4_fix进行重构
clear;clc;
warning off
close all

% ref是参考轨迹，x是射程，y是高度
global t_r u_r x_r y_r v_r gamma_r feedback w disturb
global my_nor_info
global R0 g
R0 = 10*10^3;
g = 9.81; 

my_nor_info = 1;    %使用Trajectory_normalization数据的flag

if my_nor_info
    load('Trajectory_normalization.mat')
    ref_t = Trajectory_normalization(:,1);                %控制时间序列
    ref_x = Trajectory_normalization(:,2);
    ref_y = Trajectory_normalization(:,3);
    ref_v = Trajectory_normalization(:,4);
    ref_gamma = Trajectory_normalization(:,5);
    ref_alpha = Trajectory_normalization(:,6);                %控制量u
else
    %Noop
    load('solution.mat');
end
t_r = ref_t;
x_r = ref_x;
y_r = ref_y;
v_r = ref_v;
gamma_r = ref_gamma;
u_r = ref_alpha;

feedback = 1;   %是否加控制的总开关，学习这样的写法
disturb = 1;    %是否加干扰的总开关

% 控制Maxstep和RelTol来保证较少的计算点，提高计算速度
options = odeset('Maxstep', 0.1, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %返回了这些点上的值
[T,X] = ode45(@B4yyds_Dynamics,tspan,[x_r(1),y_r(1),v_r(1),gamma_r(1)],options);
%%
t_fig = tspan*sqrt(R0/g);
y_fig = X(:,1)*R0/1000;     %km
h_fig = X(:,2)*R0/1000;     %km
v_fig = X(:,3)*sqrt(R0*g);
gama_fig = rad2deg(X(:,4));
% x y v gama
t_std_fig = t_r*sqrt(R0/g);
y_std_fig = x_r*R0/1000;     %km
h_std_fig = y_r*R0/1000;     %km
v_std_fig = v_r*sqrt(R0*g);     %m/s
gama_std_fig = rad2deg(gamma_r);    %度
%%
figure(1)
plot(y_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(y_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('射程y/km');
ylabel('高度h/km');
axis([-1 35 0 21]);
legend('标称轨迹','仿真轨迹');
title('轨迹仿真效果');

figure(2)
subplot(221)
plot(t_std_fig,y_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,y_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('射程y/km');
legend('标称轨迹','仿真轨迹');
title('轨迹仿真效果');

subplot(222)
plot(t_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('高度h/km');
legend('标称轨迹','仿真轨迹');
title('轨迹仿真效果');

subplot(223)
plot(t_std_fig,v_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,v_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('速度v/(m/s)');
legend('标称轨迹','仿真轨迹');
title('轨迹仿真效果');

subplot(224)
plot(t_std_fig,gama_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,gama_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('时间t/s');
ylabel('弹道倾角\gamma/度');
legend('标称轨迹','仿真轨迹');
title('轨迹仿真效果');
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
% %% RGB的Color不错
% figure(1)
% plot(X(:,1)*10,X(:,2)*10,'LineWidth',1)
% hold on
% plot(ref_x*10,ref_y*10,'--','Color',[0.89,0.09,0.05],'LineWidth',1.5)
% xlabel('relative lateral distance/(km)')
% ylabel('height/(km)')
% ylim([0 20.5])
% legend('CCM tracking','nominal')
% %% 看看误差
% y_fig = X(:,1)*R0/1000;
% h_fig = X(:,2)*R0/1000;
% v_fig = X(:,3)*sqrt(R0*g);
% gama_fig = rad2deg(X(:,4));