%% 预
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
%% 复合制导律
feedback = 1;   %是否加控制的总开关
disturb = 1;    %是否加干扰的总开关

% 控制Maxstep和RelTol来保证较少的计算点，提高计算速度
% step = 0.005/sqrt(R0/g);
step = 0.1;
options = odeset('Maxstep', step, 'RelTol', 1e-3);

% %24秒时制导切换为比例导引PN
% t_end = 24;
% t_end1 = t_end/sqrt(R0/g);

tspan = linspace(t_r(1) ,t_r(end) ,50); %返回了这些点上的值
[T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);

% t是行，y h v gama是列
t_fig = tspan*sqrt(R0/g);
y_fig = X(:,1)*R0/1000;     %km
h_fig = X(:,2)*R0/1000;     %km
v_fig = X(:,3)*sqrt(R0*g);
gama_fig = rad2deg(X(:,4));

% X_PN = X(50,:)';
% h_now = h_fig(end);
% alpha = interp1(t_r ,u_r, t_end1,'spline');
% 
% delta_t = 0.01/sqrt(R0/g);
% 
% for t = t_end1 : delta_t : 30/sqrt(R0/g)
%     
%     y = X_PN(1);
%     h = X_PN(2);
%     v = X_PN(3);
%     gama = X_PN(4);
%     
%     gama_wish = -atan(h/y);
%     k = -0.1*(gama-gama_wish);
%     alpha = alpha + k;
%     
%     %dynamics f，动力学
%     S   = 0.5026;                       %参考面积
%     rou = 1.225 * exp(-h*R0/7110);      %密度rou，可以不用多项式逼近
%     q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %动压
%     qf  = 0.5 * rou * v*(R0*g);               %q fake，伪动压，已约减速度v，bug标注
%     M   = v*sqrt(R0*g) / 340;                     %马赫数
%     m   = 1600;                                   %质量   
%     CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
%     L_nor = q*CL*S / (m*g);                        %升力
%     Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v
%     Cd  = 0.3042 + 0.02988*CL^2;
%     D_nor = q*Cd*S / (m*g);                        %阻力
%     
%     if disturb
%         D1_dtb = D_nor * 0.01;
%         Lf1_dtb = Lf_nor * 0.01;
%     else
%         D1_dtb = 0;
%         Lf1_dtb = 0;
%     end
%     
%     dX_PN = [-v * cos(gama);
%         v * sin(gama);
%         -D_nor - sin(gama) + D1_dtb*5;
%         Lf_nor - cos(gama)/v + Lf1_dtb*5];
%     
%     X_PN = X_PN + dX_PN * delta_t;
%     
%     if h_now < -0.01
%         break
%     end
%     
% end
%% 标称轨迹
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
%% 干扰轨迹
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
%% plot
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

%% table
Bias = {'射程偏差km';'高度偏差km';'落速m/s';'落角/度'};
Norminal = [y_std_fig(end);h_std_fig(end);v_std_fig(end);gama_std_fig(end)];
Guidance = [y_fig(end);h_fig(end);v_fig(end);gama_fig(end)];
Disturb = [y_dtb_fig(end);h_dtb_fig(end);v_dtb_fig(end);gama_dtb_fig(end)];
T = table(Bias,Norminal,Guidance,Disturb)
