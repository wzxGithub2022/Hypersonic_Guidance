%%
%-------------------------------------------------------------------------%
% normalization�������й�һ����Ϣ��t,y,h,v,gama
% M�����ù�һ������ѧ������ģ������������Ӧ�����Ƕ���ѧ����
% û��ǿ��֧�����һ�������CCM����ֱ��ʹ����ԭʼ����ѧ�ģ�����ȷʵ����ͨ�������ٽ��ܵĹ��̵õ�ʵ�ʿ�����
%-------------------------------------------------------------------------%
clc;clear;close all;
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor
global t0 t0_disp %Ϊ�˹۲�������е�ʲô�ز�

R0 = 10*10^3;
g = 9.81; 

load('Trajectory_normalization.mat')
t_nor = Trajectory_normalization(:,1);                %����ʱ������
y_nor = Trajectory_normalization(:,2);
h_nor = Trajectory_normalization(:,3);
v_nor = Trajectory_normalization(:,4);
gama_nor = Trajectory_normalization(:,5);
alpha_nor = Trajectory_normalization(:,6);                %������u

% ��һ���ĳ�ʼ״̬
y1_initial = Trajectory_normalization(1,2);
h1_initial = Trajectory_normalization(1,3);
v1_initial = Trajectory_normalization(1,4);
gama1_initial = Trajectory_normalization(1,5);
alpha1_initial = Trajectory_normalization(1,6);
state1_initial = [y1_initial, h1_initial, v1_initial, gama1_initial];         %��ֵ����35km �߶�20km �ٶ� �����ǣ��������)

%ʱ�����䣬������dt��dt1�����⻹��Ҫ���
%�������t_nor�ڷ���ģ�ΪʲôҪ�ڳ����в��ֵػָ���ʵʱ�䣬�Ѿ�����һ��������
% tspan = [ 0,  t_nor(end) * (sqrt(R0/g) ];
tspan = [ 0,  t_nor(end) ];

global step
step = 0.005;                       %����

options = odeset('RelTol',1e-6,'MaxStep',step);
[t,state1] = ode45(@BJ_Dynamitics_2D_function,tspan,state1_initial,options);  %��ά����ѧд��������

%%
%Ϊ�˻�ͼ������һ��
t_fig = t*sqrt(R0/g);       %��λ��s
y_fig = state1(:,1)/1000*R0;     %��λ��km
h_fig = state1(:,2)/1000*R0;     %��λ��km
v_fig = state1(:,3)*sqrt(R0*g);  %��λ��m/s
gama_fig = rad2deg(state1(:,4)); %��λ����

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
%Ϊ�˻�ͼ������һ��
t_norf = t_nor*sqrt(R0/g);       %��λ��s
y_norf = y_nor/1000*R0;     %��λ��km
h_norf = h_nor/1000*R0;     %��λ��km
v_norf = v_nor*sqrt(R0*g);  %��λ��m/s
gama_norf = rad2deg(gama_nor); %��λ����

figure(3),
subplot(221),
plot(t_norf,y_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,y_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');ylabel('���y/km');
legend('��ƹ켣','�ƽ�����ѧ����켣');
title('�ƽ�����ѧЧ��ͼ');
grid on;

subplot(222),
plot(t_norf,h_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');ylabel('�߶�h/km');
legend('��ƹ켣','�ƽ�����ѧ����켣');
title('�ƽ�����ѧЧ��ͼ');
grid on;

subplot(223),
plot(t_norf,v_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,v_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');ylabel('�ٶ�v/(m/s)');
legend('��ƹ켣','�ƽ�����ѧ����켣');
title('�ƽ�����ѧЧ��ͼ');
grid on;

subplot(224),
plot(t_norf,gama_norf,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,gama_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');ylabel('�������/��');
legend('��ƹ켣','�ƽ�����ѧ����켣');
title('�ƽ�����ѧЧ��ͼ');
grid on;

%%
disp('���ƫ��km');
disp(y_fig(end));
disp('�߶�ƫ��km');
disp(h_fig(end));
disp('����m/s');
disp(v_fig(end));
disp('���/��');
disp(gama_fig(end));
%%
disp('���ƫ��km');
disp(y_norf(end));
disp('�߶�ƫ��km');
disp(h_norf(end));
disp('����m/s');
disp(v_norf(end));
disp('���/��');
disp(gama_norf(end));
%% �ñƽ�����ѧ�������ӡ�Ȳ��ʵ���
% ���ƫ��km
%     4.4208
% �߶�ƫ��km
%    -1.4095
%% ��ƨ���⾫�Ȼ����Խ��ܣ���Ȼ������dot�������ܵ�4%�������������������һЩ
% ���������㣬���ڶ��Ƶ��͵����Ĳ���е��˰�
% ���ƫ��km
%     0.1386
% �߶�ƫ��km
%    -0.0220