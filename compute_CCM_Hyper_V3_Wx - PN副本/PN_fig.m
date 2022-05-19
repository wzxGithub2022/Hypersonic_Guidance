%% ����
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

feedback = 1;   %�Ƿ�ӿ��Ƶ��ܿ���
disturb = 1;    %�Ƿ�Ӹ��ŵ��ܿ���

load('X_24.mat')

%24��ʱ�Ƶ��л�Ϊ��������PN
t_end = 24;
t_end1 = t_end/sqrt(R0/g);

%%
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
    k = -0.001*(gama-gama_wish);
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
%%
figure(1)
plot(y_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
hold on
plot(y_PN_fig,h_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('���y/km');
ylabel('�߶�h/km');
axis([-1 35 0 21]);
legend('�Ƶ��켣');
title('�켣����Ч��');
grid on

figure(2)
subplot(221)
plot(t_fig,y_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
hold on
plot(t_PN_fig,y_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('���y/km');
legend('�Ƶ��켣');
title('�켣����Ч��');
grid on

subplot(222)
plot(t_fig,h_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
hold on
plot(t_PN_fig,h_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�߶�h/km');
legend('�Ƶ��켣');
title('�켣����Ч��');
grid on

subplot(223)
plot(t_fig,v_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
hold on
plot(t_PN_fig,v_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�ٶ�v/(m/s)');
legend('�Ƶ��켣');
title('�켣����Ч��');
grid on

subplot(224)
plot(t_fig,gama_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
hold on
plot(t_PN_fig,gama_PN_fig,'-','Color',[0 0.447 0.741],'LineWidth',1);
xlabel('ʱ��t/s');
ylabel('�������\gamma/��');
legend('�Ƶ��켣');
title('�켣����Ч��');
grid on

