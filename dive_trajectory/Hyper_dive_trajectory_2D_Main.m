%---------------------------------------------------%
%  最优轨迹：归一化动力学方程后二维控制落角落速问题   %
%---------------------------------------------------%
% The problem solved here is given as follows:      %
%   Minimize  alpha^2                               %
% 正负攻角都考虑到了                                 %
% 你的代价函数和升力系数阻力系数都是拿control写的     %
% cost也引入了一部分落速偏差                         %
% subject to the dynamic constraints                %
%    dy/dt =  -v1*cos(gama)                         %y是射程
%    dz/dt =  v1*sin(gama)                          %z是高度，h
%    dv/dt =  -D1/m1-sin(gama)                      %阻力致加速度，用重力加速度归一化了
%    dgama/dt=L1/(m1v1)-cos(gama)/v1                %弹道倾角
% and the boundary conditions                       %
%    y1(0) = 35000/R0                               %R0是10km，约成3.5了
%    z1(0) = 20000/R0                               %
%    v1(0) = 1750/sqrt(R0*g0)                       %sqrt(R0*g0)是313m/s
%    gama1(0)=deg2rad（-5）                         %
%    y1(t_f) = 0/R0                                 %
%    z1(t_f) = 0/R0                                 %
%    v1（t_f）=impact_velocity/sqrt(R0*g0)          %落速
%    gama1(t_f)=impact_angle                        %落角
%---------------------------------------------------%

clear; close all; clc;
%%
global R0 g0 impact_v            %落速写成global因为终端代价里用到
impact_v = 1080;                 %攻击速度 2.8Ma左右
R0 = 10*10^3;                    %R0单位：m
g0 = 9.81;
impact_angle=deg2rad(-65);       %攻击角度

%辅助数据auxdata
auxdata.g0 = 9.81; 
auxdata.S = 0.5026;              %参考面积
auxdata.R0 = 10*10^3;

%时间归一化
t10 = 0/sqrt(R0/g0); 
tf1min = 0/sqrt(R0/g0); 
tf1max = 35/sqrt(R0/g0);         %终端时间最长给到35s，sqrt(R0/g0)是31.9s

%速度与位置归一化
y10 = 35000/R0; 
z10 = 20000/R0; 
v10 = 1750/sqrt(R0*g0);  
gama0 = deg2rad(-5);
y1f = 0/R0; 
z1f = 0/R0;    
gamaf = impact_angle;              %落角严格约束
v1fmin = (impact_v-10)/sqrt(R0*g0);
v1fmax = (impact_v+10)/sqrt(R0*g0);             %落速约束+10m/s to -10m/s

%过程限制
y1min = 0/R0; y1max = 35000/R0;
z1min = 0/R0; z1max = 20000/R0;
v1min = 800/sqrt(R0*g0); v1max = 1800/sqrt(R0*g0);      %速度约束
gamamin = deg2rad(-70);  gamamax = deg2rad(25);         %弹道倾角约束
umin = deg2rad(-11);     umax = deg2rad(11);            %控制量攻角

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
iphase = 1;
bounds.phase.initialtime.lower = t10; 
bounds.phase.initialtime.upper = t10;
bounds.phase.finaltime.lower = tf1min; 
bounds.phase.finaltime.upper = tf1max;
bounds.phase.initialstate.lower = [y10,z10,v10,gama0]; 
bounds.phase.initialstate.upper = [y10,z10,v10,gama0]; 
bounds.phase.state.lower = [y1min,z1min,v1min,gamamin]; 
bounds.phase.state.upper = [y1max,z1max,v1max,gamamax]; 
bounds.phase.finalstate.lower = [y1f,z1f,v1fmin,gamaf]; 
bounds.phase.finalstate.upper = [y1f,z1f,v1fmax,gamaf]; 
bounds.phase.integral.lower = [0];          %过程代价
bounds.phase.integral.upper = [0.3];        %过程代价
bounds.phase.control.lower = umin; 
bounds.phase.control.upper = umax;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [t10; tf1max]; 
guess.phase.state   = [[y10; y1f],[z10; z1f],[v1min; v1max],[gamamin;gamamax]];
guess.phase.control = [-deg2rad(13); deg2rad(13)];      %过程限制11，这里猜测给到13，稍微放宽
guess.phase.integral = 0.1;                             %过程代价猜测0.1

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%
%-------------------------------------------------------------------------%
setup.name = 'Hyper_dive_trajectory_2D_Problem';
setup.functions.continuous = @Hyper_dive_trajectory_2D_Continuous;
setup.functions.endpoint = @Hyper_dive_trajectory_2D_Endpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'snopt';                     %已经集成好的NPL解算器
setup.derivatives.supplier = 'sparseCD';        %sparseFD sparseBD or sparseCD，default sparseFD
setup.derivatives.derivativelevel = 'second';   %first or second，default first
setup.mesh.method = 'hp1';                      %hp or hp1,default hp1
%setup.mesh.tolerance = 1e-4;
setup.mesh.tolerance = 1e-6;
setup.mesh.maxiteration = 40;   %最大迭代次数
setup.mesh.colpointsmin = 4;
setup.mesh.colpointsmax = 10;
setup.mesh.phase.colpoints = 4*ones(1,10);
setup.mesh.phase.fraction =  0.1*ones(1,10);
setup.method = 'RPMintegration';

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOP2 ---------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
solution = output.result.solution;

%--------------------------------------------------------------------------%
%------------------------------- Plot Solution ----------------------------%
%--------------------------------------------------------------------------%
%注意：时间序列步长不是均匀的
%最终运行时间,这个end很灵性诶，不打分号是为了在命令行显示出来，下同
t_final = solution.phase(1).time(end)*sqrt(R0/g0)       

%误差统计
final_impact_angle = solution.phase(1).state(end,4);              %最终落角rad
error_angle = rad2deg(final_impact_angle-impact_angle)            %落角误差(单位：度)
final_impact_v = solution.phase(1).state(end,3)*sqrt(R0*g0);      %最终落速
error_velocity = final_impact_v-impact_v                          %落速误差
miss_distance = sqrt((solution.phase(1).state(end,2)*R0)^2+(solution.phase(1).state(end,1)*R0)^2) %脱靶量，横程+纵程

%% 随附绘图
%为了绘图，反归一化
t_fig = solution.phase(1).time*sqrt(R0/g0);       %单位：s
y_fig = solution.phase(1).state(:,1)/1000*R0;     %单位：km
h_fig = solution.phase(1).state(:,2)/1000*R0;     %单位：km
v_fig = solution.phase(1).state(:,3)*sqrt(R0*g0); %单位：m/s
gama_fig = rad2deg(solution.phase(1).state(:,4)); %单位：度
u_fig = rad2deg(solution.phase(1).control);       %单位：度

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
plot(t_fig,u_fig,'LineWidth',2);
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

%% 动力学仿真验证
% t_seq = solution.phase(1).time*sqrt(R0/g0);       %单位：s
% u_seq = solution.phase(1).control;                %单位：弧度
% time_control_sequence = [t_seq,u_seq];
% save('Trajectory_control_information.mat','time_control_sequence');     %轨迹仿真的控制信息
% 
% y_fig = solution.phase(1).state(:,1)/1000*R0;     %单位：km
% h_fig = solution.phase(1).state(:,2)/1000*R0;     %单位：km
% trajectory_location_sequence = [y_fig,h_fig];
% save('Trajectory_opt_information.mat','trajectory_location_sequence');  %轨迹仿真的位置信息


%% 线性化时计算需求
% % 半可视化数据
% % t_fig 单位：s
% % y_fig 单位：km
% % h_fig 单位：km
% % v_fig 单位：m/s
% gama_seq = solution.phase(1).state(:,4);    %单位：弧度
% % u_seq 单位：弧度
% Trajectory_all_information = [t_fig,y_fig,h_fig,v_fig,gama_seq,u_seq];
% save('Trajectory_all_information.mat','Trajectory_all_information');

% normalization，归一化数据
t_nor = solution.phase(1).time;
y_nor = solution.phase(1).state(:,1);     
h_nor = solution.phase(1).state(:,2);     
v_nor = solution.phase(1).state(:,3);     
gama_nor = solution.phase(1).state(:,4); 
alpha_nor = solution.phase(1).control;       
Trajectory_normalization = [t_nor, y_nor, h_nor, v_nor, gama_nor, alpha_nor];
save('Trajectory_normalization.mat','Trajectory_normalization');