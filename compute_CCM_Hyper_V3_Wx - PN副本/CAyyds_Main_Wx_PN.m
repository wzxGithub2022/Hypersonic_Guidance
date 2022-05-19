% %% ����
% %���������㷨��ܣ���CA�����ع�
% %�ģ������Եĵ�����ȥ��msspoly�����Ҳ���������CCM�����ǵ��޶���һ��С�����ڣ������Ի�Ҳ��𲻴�
% clear;clc;
% warning off
% close all
% 
% % ref�ǲο��켣��x����̣�y�Ǹ߶�
% global t_r u_r y_r h_r v_r gamma_r feedback disturb
% global R0 g
% R0 = 10*10^3;
% g = 9.81; 
% 
% load('Trajectory_normalization.mat')
% t_r = Trajectory_normalization(:,1);
% y_r = Trajectory_normalization(:,2);
% h_r = Trajectory_normalization(:,3);
% v_r = Trajectory_normalization(:,4);
% gamma_r = Trajectory_normalization(:,5);
% u_r = Trajectory_normalization(:,6);
% %% Wx�Ƶ�
% feedback = 1;   %�Ƿ�ӿ��Ƶ��ܿ���
% disturb = 1;    %�Ƿ�Ӹ��ŵ��ܿ���
% 
% % ����Maxstep��RelTol����֤���ٵļ���㣬��߼����ٶ�
% % step = 0.005/sqrt(R0/g);
% step = 0.1;
% options = odeset('Maxstep', step, 'RelTol', 1e-3);
% 
% %24��ʱ�Ƶ��л�Ϊ��������PN
% t_end = 24;
% t_end1 = t_end/sqrt(R0/g);
% 
% tspan = linspace(t_r(1) ,t_end1 ,50); %��������Щ���ϵ�ֵ
% [T,X] = ode45(@CAyyds_Dynamics,tspan,[y_r(1),h_r(1),v_r(1),gamma_r(1)],options);
% 
% t_fig = tspan*sqrt(R0/g);
% y_fig = X(:,1)*R0/1000;     %km
% h_fig = X(:,2)*R0/1000;     %km
% v_fig = X(:,3)*sqrt(R0*g);
% gama_fig = rad2deg(X(:,4));
%% ���Ʒ���� ���� Wx�Ƶ� ���ɼ��
clear;clc;
warning off
close all
load('X_24.mat')
global R0 g
%% PN����
X_PN = X(50,:)';
h_now = h_fig(end);
alpha = interp1(t_r ,u_r, t_end1,'spline');

delta_t = 0.01/sqrt(R0/g);

i=1;

for t = t_end1 : delta_t : 30/sqrt(R0/g)
    
    t*sqrt(R0/g)
    
    t_PN(i) = t;
    y_PN(i) = X_PN(1);
    h_PN(i) = X_PN(2);
    v_PN(i) = X_PN(3);
    gama_PN(i) = X_PN(4);
    alpha_PN(i) = alpha;
    i=i+1;
    
    y = X_PN(1);
    h = X_PN(2);
    v = X_PN(3);
    gama = X_PN(4);
    
    gama_wish = -atan(h/y);
    k = -0.00255*(gama-gama_wish);
    
%     k = -0.001*(gama-gama_wish)*norm([y*R0,h*R0]);
%     %����
%     if norm([y*R0,h*R0]) < 1000
%         k=0;
%     end

%     Az = -atan(h/y);
%     k = 0.002 * ( Az - deg2rad(-65) );

    %����ûë�����ø�ƫ�ñ�������
    k_lim = 0.001;
    if k > k_lim
        k = k_lim;
    end
    if k < -k_lim
        k = -k_lim;
    end
    
    alpha = alpha + k;
    
    %dynamics f������ѧ
    S   = 0.5026;                       %�ο����
    rou = 1.225 * exp(-h*R0/7110);      %�ܶ�rou�����Բ��ö���ʽ�ƽ�
    q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
    qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
    M   = v*sqrt(R0*g) / 340;                     %�����
    m   = 1600;                                   %����   
    CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
    L_nor = q*CL*S / (m*g);                        %����
    Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v
    Cd  = 0.3042 + 0.02988*CL^2;
    D_nor = q*Cd*S / (m*g);                        %����
    
    if disturb
        D1_dtb = D_nor * 0.01;
        Lf1_dtb = Lf_nor * 0.01;
    else
        D1_dtb = 0;
        Lf1_dtb = 0;
    end
    
    dX_PN = [-v * cos(gama);
        v * sin(gama);
        -D_nor - sin(gama) + D1_dtb*5;
        Lf_nor - cos(gama)/v + Lf1_dtb*5];
    
    X_PN = X_PN + dX_PN * delta_t;
    
    if h*R0 < 1
        break
    end
    
end

t_PN_fig = t_PN*sqrt(R0/g);
y_PN_fig = y_PN*R0/1000;
h_PN_fig = h_PN*R0/1000;
v_PN_fig = v_PN*sqrt(R0*g);
gama_PN_fig = rad2deg(gama_PN);
alpha_PN_fig = rad2deg(alpha_PN);
%% ��ƹ켣
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
%% ���Ź켣
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
%% plot
figure(1)
plot(y_std_fig,h_std_fig,'--','Color',[0.85 0.325 0.098],'LineWidth',2);
hold on
plot(y_dtb_fig,h_dtb_fig,':','Color',[0.929 0.694 0.125],'LineWidth',1);
hold on
plot(y_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
hold on
plot(y_PN_fig,h_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
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
hold on
plot(t_PN_fig,y_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
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
hold on
plot(t_PN_fig,h_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
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
hold on
plot(t_PN_fig,v_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
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
hold on
plot(t_PN_fig,gama_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�������\gamma/��');
legend('��ƹ켣','���Ź켣','�Ƶ��켣');
title('�켣����Ч��');
grid on

%% table
Bias = {'���ƫ��km';'�߶�ƫ��km';'����m/s';'���/��'};
Norminal = [y_std_fig(end);h_std_fig(end);v_std_fig(end);gama_std_fig(end)];
Guidance = [y_PN_fig(end);h_PN_fig(end);v_PN_fig(end);gama_PN_fig(end)];
Disturb = [y_dtb_fig(end);h_dtb_fig(end);v_dtb_fig(end);gama_dtb_fig(end)];
T = table(Bias,Norminal,Guidance,Disturb)
