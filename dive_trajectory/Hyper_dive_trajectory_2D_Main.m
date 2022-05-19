%---------------------------------------------------%
%  ���Ź켣����һ������ѧ���̺��ά���������������   %
%---------------------------------------------------%
% The problem solved here is given as follows:      %
%   Minimize  alpha^2                               %
% �������Ƕ����ǵ���                                 %
% ��Ĵ��ۺ���������ϵ������ϵ��������controlд��     %
% costҲ������һ��������ƫ��                         %
% subject to the dynamic constraints                %
%    dy/dt =  -v1*cos(gama)                         %y�����
%    dz/dt =  v1*sin(gama)                          %z�Ǹ߶ȣ�h
%    dv/dt =  -D1/m1-sin(gama)                      %�����¼��ٶȣ����������ٶȹ�һ����
%    dgama/dt=L1/(m1v1)-cos(gama)/v1                %�������
% and the boundary conditions                       %
%    y1(0) = 35000/R0                               %R0��10km��Լ��3.5��
%    z1(0) = 20000/R0                               %
%    v1(0) = 1750/sqrt(R0*g0)                       %sqrt(R0*g0)��313m/s
%    gama1(0)=deg2rad��-5��                         %
%    y1(t_f) = 0/R0                                 %
%    z1(t_f) = 0/R0                                 %
%    v1��t_f��=impact_velocity/sqrt(R0*g0)          %����
%    gama1(t_f)=impact_angle                        %���
%---------------------------------------------------%

clear; close all; clc;
%%
global R0 g0 impact_v            %����д��global��Ϊ�ն˴������õ�
impact_v = 1080;                 %�����ٶ� 2.8Ma����
R0 = 10*10^3;                    %R0��λ��m
g0 = 9.81;
impact_angle=deg2rad(-65);       %�����Ƕ�

%��������auxdata
auxdata.g0 = 9.81; 
auxdata.S = 0.5026;              %�ο����
auxdata.R0 = 10*10^3;

%ʱ���һ��
t10 = 0/sqrt(R0/g0); 
tf1min = 0/sqrt(R0/g0); 
tf1max = 35/sqrt(R0/g0);         %�ն�ʱ�������35s��sqrt(R0/g0)��31.9s

%�ٶ���λ�ù�һ��
y10 = 35000/R0; 
z10 = 20000/R0; 
v10 = 1750/sqrt(R0*g0);  
gama0 = deg2rad(-5);
y1f = 0/R0; 
z1f = 0/R0;    
gamaf = impact_angle;              %����ϸ�Լ��
v1fmin = (impact_v-10)/sqrt(R0*g0);
v1fmax = (impact_v+10)/sqrt(R0*g0);             %����Լ��+10m/s to -10m/s

%��������
y1min = 0/R0; y1max = 35000/R0;
z1min = 0/R0; z1max = 20000/R0;
v1min = 800/sqrt(R0*g0); v1max = 1800/sqrt(R0*g0);      %�ٶ�Լ��
gamamin = deg2rad(-70);  gamamax = deg2rad(25);         %�������Լ��
umin = deg2rad(-11);     umax = deg2rad(11);            %����������

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
bounds.phase.integral.lower = [0];          %���̴���
bounds.phase.integral.upper = [0.3];        %���̴���
bounds.phase.control.lower = umin; 
bounds.phase.control.upper = umax;

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [t10; tf1max]; 
guess.phase.state   = [[y10; y1f],[z10; z1f],[v1min; v1max],[gamamin;gamamax]];
guess.phase.control = [-deg2rad(13); deg2rad(13)];      %��������11������²����13����΢�ſ�
guess.phase.integral = 0.1;                             %���̴��۲²�0.1

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%
%-------------------------------------------------------------------------%
setup.name = 'Hyper_dive_trajectory_2D_Problem';
setup.functions.continuous = @Hyper_dive_trajectory_2D_Continuous;
setup.functions.endpoint = @Hyper_dive_trajectory_2D_Endpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'snopt';                     %�Ѿ����ɺõ�NPL������
setup.derivatives.supplier = 'sparseCD';        %sparseFD sparseBD or sparseCD��default sparseFD
setup.derivatives.derivativelevel = 'second';   %first or second��default first
setup.mesh.method = 'hp1';                      %hp or hp1,default hp1
%setup.mesh.tolerance = 1e-4;
setup.mesh.tolerance = 1e-6;
setup.mesh.maxiteration = 40;   %����������
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
%ע�⣺ʱ�����в������Ǿ��ȵ�
%��������ʱ��,���end��������������ֺ���Ϊ������������ʾ��������ͬ
t_final = solution.phase(1).time(end)*sqrt(R0/g0)       

%���ͳ��
final_impact_angle = solution.phase(1).state(end,4);              %�������rad
error_angle = rad2deg(final_impact_angle-impact_angle)            %������(��λ����)
final_impact_v = solution.phase(1).state(end,3)*sqrt(R0*g0);      %��������
error_velocity = final_impact_v-impact_v                          %�������
miss_distance = sqrt((solution.phase(1).state(end,2)*R0)^2+(solution.phase(1).state(end,1)*R0)^2) %�Ѱ��������+�ݳ�

%% �渽��ͼ
%Ϊ�˻�ͼ������һ��
t_fig = solution.phase(1).time*sqrt(R0/g0);       %��λ��s
y_fig = solution.phase(1).state(:,1)/1000*R0;     %��λ��km
h_fig = solution.phase(1).state(:,2)/1000*R0;     %��λ��km
v_fig = solution.phase(1).state(:,3)*sqrt(R0*g0); %��λ��m/s
gama_fig = rad2deg(solution.phase(1).state(:,4)); %��λ����
u_fig = rad2deg(solution.phase(1).control);       %��λ����

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
plot(t_fig,u_fig,'LineWidth',2);
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

%% ����ѧ������֤
% t_seq = solution.phase(1).time*sqrt(R0/g0);       %��λ��s
% u_seq = solution.phase(1).control;                %��λ������
% time_control_sequence = [t_seq,u_seq];
% save('Trajectory_control_information.mat','time_control_sequence');     %�켣����Ŀ�����Ϣ
% 
% y_fig = solution.phase(1).state(:,1)/1000*R0;     %��λ��km
% h_fig = solution.phase(1).state(:,2)/1000*R0;     %��λ��km
% trajectory_location_sequence = [y_fig,h_fig];
% save('Trajectory_opt_information.mat','trajectory_location_sequence');  %�켣�����λ����Ϣ


%% ���Ի�ʱ��������
% % ����ӻ�����
% % t_fig ��λ��s
% % y_fig ��λ��km
% % h_fig ��λ��km
% % v_fig ��λ��m/s
% gama_seq = solution.phase(1).state(:,4);    %��λ������
% % u_seq ��λ������
% Trajectory_all_information = [t_fig,y_fig,h_fig,v_fig,gama_seq,u_seq];
% save('Trajectory_all_information.mat','Trajectory_all_information');

% normalization����һ������
t_nor = solution.phase(1).time;
y_nor = solution.phase(1).state(:,1);     
h_nor = solution.phase(1).state(:,2);     
v_nor = solution.phase(1).state(:,3);     
gama_nor = solution.phase(1).state(:,4); 
alpha_nor = solution.phase(1).control;       
Trajectory_normalization = [t_nor, y_nor, h_nor, v_nor, gama_nor, alpha_nor];
save('Trajectory_normalization.mat','Trajectory_normalization');