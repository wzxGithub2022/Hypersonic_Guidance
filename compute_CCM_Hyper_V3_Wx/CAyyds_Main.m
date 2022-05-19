%���������㷨��ܣ���CA�����ع�
%�ģ������Եĵ�����ȥ��msspoly�����Ҳ���������CCM�����ǵ��޶���һ��С�����ڣ������Ի�Ҳ��𲻴�
clear;clc;
warning off
close all

% ref�ǲο��켣��x����̣�y�Ǹ߶�
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
feedback = 1;   %�Ƿ�ӿ��Ƶ��ܿ���
disturb = 1;    %�Ƿ�Ӹ��ŵ��ܿ���

% ����Maxstep��RelTol����֤���ٵļ���㣬��߼����ٶ�
% step = 0.005/sqrt(R0/g);
step = 0.1;
options = odeset('Maxstep', step, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %��������Щ���ϵ�ֵ
[T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);

t_fig = tspan*sqrt(R0/g);
y_fig = X(:,1)*R0/1000;     %km
h_fig = X(:,2)*R0/1000;     %km
v_fig = X(:,3)*sqrt(R0*g);
gama_fig = rad2deg(X(:,4));
%%
feedback = 0;   %�Ƿ�ӿ��Ƶ��ܿ���
disturb = 0;    %�Ƿ�Ӹ��ŵ��ܿ���

% ����Maxstep��RelTol����֤���ٵļ���㣬��߼����ٶ�
% step = 0.005/sqrt(R0/g);
step = 0.1;
options = odeset('Maxstep', step, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %��������Щ���ϵ�ֵ
[T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);

t_std_fig = tspan*sqrt(R0/g);
y_std_fig = X(:,1)*R0/1000;     %km
h_std_fig = X(:,2)*R0/1000;     %km
v_std_fig = X(:,3)*sqrt(R0*g);
gama_std_fig = rad2deg(X(:,4));
%%
feedback = 0;   %�Ƿ�ӿ��Ƶ��ܿ���
disturb = 1;    %�Ƿ�Ӹ��ŵ��ܿ���

% ����Maxstep��RelTol����֤���ٵļ���㣬��߼����ٶ�
% step = 0.005/sqrt(R0/g);
step = 0.1;
options = odeset('Maxstep', step, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %��������Щ���ϵ�ֵ
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
xlabel('���y/km');
ylabel('�߶�h/km');
axis([-1 35 0 21]);
legend('��ƹ켣','���Ź켣','�Ƶ��켣');
title('�켣����Ч��');
grid on

figure(2)
subplot(221)
plot(t_std_fig,y_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,y_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,y_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('���y/km');
legend('��ƹ켣','���Ź켣','�Ƶ��켣');
title('�켣����Ч��');
grid on

subplot(222)
plot(t_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,h_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�߶�h/km');
legend('��ƹ켣','���Ź켣','�Ƶ��켣');
title('�켣����Ч��');
grid on

subplot(223)
plot(t_std_fig,v_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,v_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,v_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�ٶ�v/(m/s)');
legend('��ƹ켣','���Ź켣','�Ƶ��켣');
title('�켣����Ч��');
grid on

subplot(224)
plot(t_std_fig,gama_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_dtb_fig,gama_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(t_fig,gama_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�������\gamma/��');
legend('��ƹ켣','���Ź켣','�Ƶ��켣');
title('�켣����Ч��');
grid on

%%
Bias = {'���ƫ��km';'�߶�ƫ��km';'����m/s';'���/��'};
Norminal = [y_std_fig(end);h_std_fig(end);v_std_fig(end);gama_std_fig(end)];
Guidance = [y_fig(end);h_fig(end);v_fig(end);gama_fig(end)];
Disturb = [y_dtb_fig(end);h_dtb_fig(end);v_dtb_fig(end);gama_dtb_fig(end)];
T = table(Bias,Norminal,Guidance,Disturb)
%%
% disp('�Ƶ�-���ƫ��km');
% disp(y_fig(end));
% disp('�Ƶ�-�߶�ƫ��km');
% disp(h_fig(end));
% disp('�Ƶ�-����m/s');
% disp(v_fig(end));
% disp('�Ƶ�-���/��');
% disp(gama_fig(end));
% 
% disp('�Ŷ�-���ƫ��km');
% disp(y_dtb_fig(end));
% disp('�Ŷ�-�߶�ƫ��km');
% disp(h_dtb_fig(end));
% disp('�Ŷ�-����m/s');
% disp(v_dtb_fig(end));
% disp('�Ŷ�-���/��');
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
% % title('��ά�켣����')
% %%
% % ���ƫ��km
% %    -0.0037
% % �߶�ƫ��km
% %    -0.0475
% % ����m/s
% %    1.0801e+03
% % ���/��
% %   -65.4138