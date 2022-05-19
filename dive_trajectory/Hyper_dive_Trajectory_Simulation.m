%-------------------------------------------------------------------------%
%------------------ plot the trajectory with ode -------------------------%
%------------------ 2D_plain(y,z)-----------------------------------------%
%------------------ target location  (0,0) -------------------------------%
%------------------ initial location (y0,z0) -----------------------------%
%-------------------------------------------------------------------------%
clc;clear;close all;
%%
global t_u u

R0 = 10*10^3;
g0 = 9.81; 

load('Trajectory_control_information.mat')
t_u = time_control_sequence(:,1);                %����ʱ������
u   = time_control_sequence(:,2);                %������u

impact_angle = deg2rad(-65);                     %���Լ��
tspan = [0 t_u(end)];                            %ʱ������
y0 = [35*10^3 20*10^3 1750 deg2rad(-5)];         %��ֵ����35km �߶�20km �ٶ� �����ǣ��������)
step = 0.005;                                    %����
options = odeset('RelTol',1e-6,'MaxStep',step);
[t,y] = ode45(@Hyper_dive_Dynamitics_2D_function,tspan,y0,options);  %��ά����ѧд��������

%%
%��ͼ
y_sim = y(:,1)/1000;
h_sim = y(:,2)/1000;

load('Trajectory_opt_information.mat')
y_opt = trajectory_location_sequence(:,1);
h_opt = trajectory_location_sequence(:,2);

figure(1)
hold on,plot(y_opt,h_opt,'-*','Color',[0.85 0.325 0.098],'LineWidth',1);
hold on,plot(y_sim,h_sim,'--','Color',[0 0.447 0.741],'LineWidth',2);
xlabel('����y/km');ylabel('�߶�h/km');
legend('�Ż�����켣','����ѧ����켣');
title('����ѧ�������Ż�����켣�Ա�');
grid on;

%%
%���ͳ��
fprintf('�ն˺��̣�%d (m)\n',  y(end,1));
fprintf('�ն˸߶ȣ�%d (m)\n',  y(end,2));
fprintf('�ն����٣�%d (m/s)\n',y(end,3));
fprintf('����ƫ�%d (m/s)\n',y(end,3)-1080);
fprintf('�ն���ǣ�%d (��)\n',rad2deg(y(end,4)));
fprintf('���ƫ�%d (��)\n',rad2deg(y(end,4))-(-65));
%%
% �ն˺��̣�-7.684552e-02 (m)
% �ն˸߶ȣ�3.606829e-02 (m)
% �ն����٣�1.076287e+03 (m/s)
% �ն���ǣ�-6.499940e+01 (��)