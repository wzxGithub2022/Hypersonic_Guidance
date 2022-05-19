%%
%-------------------------------------------------------------------------%
% normalization�������й�һ����Ϣ��t,y,h,v,gama
% M�����ù�һ������ѧ������ģ������������Ӧ�����Ƕ���ѧ����
% û��ǿ��֧�����һ�������CCM����ֱ��ʹ����ԭʼ����ѧ�ģ�����ȷʵ����ͨ�������ٽ��ܵĹ��̵õ�ʵ�ʿ�����
% R0��10km��sqrt(R0*g0)��313m/s��sqrt(R0/g0)��31.9s
%-------------------------------------------------------------------------%
clc;clear;close all;
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor B_nor df_nor 
global t0 t0_disp %Ϊ�˹۲�������е�ʲô�ز�
global feedback_QP
global A_fig b_fig
A_fig = [];
b_fig = [];

feedback_QP = 1;

R0 = 10*10^3;
g = 9.81; 

load('Trajectory_normalization.mat')
load('CCM_upper_B4.mat')
t_nor = Trajectory_normalization(:,1);                %����ʱ������
y_nor = Trajectory_normalization(:,2);
h_nor = Trajectory_normalization(:,3);
v_nor = Trajectory_normalization(:,4);
gama_nor = Trajectory_normalization(:,5);
alpha_nor = Trajectory_normalization(:,6);                %������u
CCM_nor = CCM_upper;

load('CCM_Dynamitics_B.mat');
load('CCM_Dynamitics_df.mat');
B_nor = B_mat_value;
df_nor = df_mat_value;

% ��һ���ĳ�ʼ״̬
y1_initial = Trajectory_normalization(1,2);
h1_initial = Trajectory_normalization(1,3);
v1_initial = Trajectory_normalization(1,4);
gama1_initial = Trajectory_normalization(1,5);
alpha1_initial = Trajectory_normalization(1,6);
state1_initial = [y1_initial, h1_initial, v1_initial, gama1_initial, alpha1_initial];         %��ֵ����35km �߶�20km �ٶ� �����ǣ��������)

tspan = [0 t_nor(end)];                            %ʱ������
global step
step = 0.005/(sqrt(R0/g)) * 100;                       %�����������Ŵ���ܷ���Ҳ����ˣ��ܺ�

options = odeset('RelTol',1e-3,'MaxStep',step);
[t,state1] = ode45(@Hyper_dive_Dynamitics_2D_function,tspan,state1_initial,options);  %��ά����ѧд��������

%%
%Ϊ�˻�ͼ������һ��
t_fig = t*sqrt(R0/g);       %��λ��s
y_fig = state1(:,1)/1000*R0;     %��λ��km
h_fig = state1(:,2)/1000*R0;     %��λ��km
v_fig = state1(:,3)*sqrt(R0*g);  %��λ��m/s
gama_fig = rad2deg(state1(:,4)); %��λ����
alpha_fig = rad2deg(state1(:,5)); %��λ����

figure(1),
plot(t_fig,h_fig,'-ob',t_fig,y_fig,'-*r');
xlabel('ʱ��t/s');ylabel('�߶������/km');legend('�߶�h','���y');
title('�߶��������ʱ��ı仯����');
grid on

figure(2),subplot(221),
plot(y_fig,h_fig,'linewidth',2);
xlabel('���y/km');ylabel('�߶�h/km');
title('��άƽ�����Ź켣');
grid on

subplot(222),
plot(t_fig,alpha_fig,'LineWidth',2);
xlabel('ʱ��t/s');ylabel('����\alpha/��');
title('������ʱ��ı仯����');
grid on

subplot(223),
plot(t_fig,v_fig,'linewidth',3);
xlabel('ʱ��t/s');ylabel('�ٶ�v/(m/s)');
title('�ٶ���ʱ��ı仯����');
grid on

subplot(224),
plot(t_fig,gama_fig,'linewidth',3);
xlabel('ʱ��t/s');ylabel('�������/��');
title('���������ʱ��ı仯����');
grid on

%%
disp('���ƫ��km');
disp(y_fig(end));
disp('�߶�ƫ��km');
disp(h_fig(end));

%%
figure(3);
subplot(121),plot(A_fig);title('A_fig');
subplot(122),plot(b_fig);title('b_fig');
