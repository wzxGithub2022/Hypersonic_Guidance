%Ϊ����߳���Ŀɱ��Ժ�����Ч�ʣ����������㷨��ܣ���B4_fix�����ع�
clear;clc;
warning off
close all

% ref�ǲο��켣��x����̣�y�Ǹ߶�
global t_r u_r x_r y_r v_r gamma_r feedback w disturb
global my_nor_info
global R0 g
R0 = 10*10^3;
g = 9.81; 

my_nor_info = 1;    %ʹ��Trajectory_normalization���ݵ�flag

if my_nor_info
    load('Trajectory_normalization.mat')
    ref_t = Trajectory_normalization(:,1);                %����ʱ������
    ref_x = Trajectory_normalization(:,2);
    ref_y = Trajectory_normalization(:,3);
    ref_v = Trajectory_normalization(:,4);
    ref_gamma = Trajectory_normalization(:,5);
    ref_alpha = Trajectory_normalization(:,6);                %������u
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

feedback = 1;   %�Ƿ�ӿ��Ƶ��ܿ��أ�ѧϰ������д��
disturb = 1;    %�Ƿ�Ӹ��ŵ��ܿ���

% ����Maxstep��RelTol����֤���ٵļ���㣬��߼����ٶ�
options = odeset('Maxstep', 0.1, 'RelTol', 1e-3);

tspan = linspace(t_r(1) ,t_r(end) ,50); %��������Щ���ϵ�ֵ
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
gama_std_fig = rad2deg(gamma_r);    %��
%%
figure(1)
plot(y_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(y_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('���y/km');
ylabel('�߶�h/km');
axis([-1 35 0 21]);
legend('��ƹ켣','����켣');
title('�켣����Ч��');

figure(2)
subplot(221)
plot(t_std_fig,y_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,y_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('���y/km');
legend('��ƹ켣','����켣');
title('�켣����Ч��');

subplot(222)
plot(t_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�߶�h/km');
legend('��ƹ켣','����켣');
title('�켣����Ч��');

subplot(223)
plot(t_std_fig,v_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,v_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�ٶ�v/(m/s)');
legend('��ƹ켣','����켣');
title('�켣����Ч��');

subplot(224)
plot(t_std_fig,gama_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(t_fig,gama_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�������\gamma/��');
legend('��ƹ켣','����켣');
title('�켣����Ч��');
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
% %% RGB��Color����
% figure(1)
% plot(X(:,1)*10,X(:,2)*10,'LineWidth',1)
% hold on
% plot(ref_x*10,ref_y*10,'--','Color',[0.89,0.09,0.05],'LineWidth',1.5)
% xlabel('relative lateral distance/(km)')
% ylabel('height/(km)')
% ylim([0 20.5])
% legend('CCM tracking','nominal')
% %% �������
% y_fig = X(:,1)*R0/1000;
% h_fig = X(:,2)*R0/1000;
% v_fig = X(:,3)*sqrt(R0*g);
% gama_fig = rad2deg(X(:,4));