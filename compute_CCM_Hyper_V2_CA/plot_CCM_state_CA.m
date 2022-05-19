%%
load('state_CCM_CA.mat');
x = linspace(1,46,46);
%%
% alpha_w = max(sigma_ThBw(:));       %(:)��ʾ��������Ԫ�ز����б�ʾ
% d_bar = alpha_w/lambda;             %��d_barû�г�w�Ͻ磬������J_CCMҲ�����ٳ�w�Ͻ��ˣ�һ����λ
% disp('d_bar'); 
% disp(d_bar);
figure(1)
plot(x,state_CCM(1,:));
title('d bar');
%%
% disp('Control:'); 
% disp(max(d_bar*delta_u(:)));        %���ǿ����Ͻ�
figure(2)
plot(x,state_CCM(2,:),'-o','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('ѡȡ�ĵ�i��״̬��');
ylabel('�����Ͻ�');
title('����ѡȡ���µļ�����������');
grid on
%%
% disp('W:'); 
% disp(min(min(min(eig_W(:,:,:,1)))));      %С��С
% disp(max(max(max(eig_W(:,:,:,2)))));      %���д�
figure(3)
subplot(121),plot(x,state_CCM(3,:));
title('W min');
subplot(122),plot(x,state_CCM(4,:));
title('W max');
%%
% disp('min eig CCM:'); 
% disp(min(eig_CCM_min(:)));
% disp('max eig CCM:'); 
% disp(max(eig_CCM_max(:)));
figure(4)
subplot(121),plot(x,state_CCM(5,:));
title('R-CCM eig min');
subplot(122),plot(x,state_CCM(6,:));
title('R-CCM eig max');
%%
% disp('euc_bounds');
% disp(d_bar*sqrt(diag(W_upper)));
figure(5)
subplot(231),plot(x,state_CCM(7,:));
title('euc bounds 1');
subplot(232),plot(x,state_CCM(8,:));
title('euc bounds 2');
subplot(233),plot(x,state_CCM(9,:));
title('euc bounds 3');
subplot(234),plot(x,state_CCM(10,:));
title('euc bounds 4');
%%
load('Trajectory_normalization.mat')

global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

t = Trajectory_normalization(:,1)*sqrt(R0/g);

figure(6)
plot(t,state_CCM(1,:),'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('ʱ��t/s');
ylabel('Tube \alpha_w/\lambda');
title('���ڵ�λ�Ŷ���Tube');
grid on

figure(7)
subplot(221),
plot(t,state_CCM(7,:)*R0,'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('ʱ��t/s');
ylabel('Tube ���y/m');
title('���ڵ�λ�Ŷ���Tube');
grid on

subplot(222),
plot(t,state_CCM(8,:)*R0,'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('ʱ��t/s');
ylabel('Tube �߶�h/m');
title('���ڵ�λ�Ŷ���Tube');
grid on

subplot(223),
plot(t,state_CCM(9,:)*sqrt(R0*g),'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('ʱ��t/s');
ylabel('Tube �ٶ�v/(m/s)');
title('���ڵ�λ�Ŷ���Tube');
grid on

subplot(224),
plot(t,rad2deg( state_CCM(10,:) ),'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('ʱ��t/s');
ylabel('Tube �������\gamma/��');
title('���ڵ�λ�Ŷ���Tube');
grid on
%%
figure(8)
plot(t,( state_CCM(7,:) - state_CCM(8,:) )*R0,'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('ʱ��t/s');
ylabel('�����߶ȵ�Tubeƫ��/m');
title('���ڵ�λ�Ŷ���Tube');
grid on
%%
clc
x = Trajectory_normalization(:,2)*R0/1000;
y = zeros(length(x));
z = Trajectory_normalization(:,3)*R0/1000;
figure(9)
L0 = plot3(x,y,z,'--*','Color',[0 0.447 0.741],'Linewidth',1);

% legend(L0 ,'��ƹ켣')
num = length(x);

for i=1:1
    x0 = x(i);
    y0 = y(i);
    z0 = z(i);
    r = state_CCM(7,i)*R0/1000;
%     theta = linspace(-pi,pi,length(x));
    theta = linspace(-pi,pi,100);
    xr = x0 + zeros(length(theta));
    yr = y0 + r*cos(theta);
    zr = z0 + r*sin(theta);
    hold on
    L1 = plot3(xr,yr,zr,'-','Color',[0.85 0.325 0.098],'Linewidth',2);
end
% ah=axes('position',get(gca,'position'),'visible','off');
% legend(L1 ,'Tube')
for i=2:num
    x0 = x(i);
    y0 = y(i);
    z0 = z(i);
    r = state_CCM(7,i)*R0/1000;
    theta = linspace(-pi,pi,100);
    xr = x0 + zeros(length(theta));
    yr = y0 + r*cos(theta);
    zr = z0 + r*sin(theta);
    hold on
    plot3(xr,yr,zr,'-','Color',[0.85 0.325 0.098],'Linewidth',2);
end
hold off

axis([0 35 -10 10 -2 20.5]);

xlabel('��̷���/km');
ylabel('��չ����/km');
zlabel('�߶ȷ���/km');
title('���ڱ�ƹ켣����ά��չTube');
grid on

% legend(L0 ,'��ƹ켣','Location','northeast')
% legend('boxoff')
% ah=axes('position',get(gca,'position'),'visible','off');
% legend(ah,L1 ,'Tube','Location','northeast')
% legend('boxoff')

% legend([L0 L1],{'��ƹ켣','Tube'})
%%
clear;clc;close all

load('Trajectory_normalization.mat')

global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

x = Trajectory_normalization(:,2)*R0/1000;
y = zeros(length(x));
z = Trajectory_normalization(:,3)*R0/1000;
figure(10)

plot3(x,y,z,'--*','Color',[0 0.447 0.741],'Linewidth',1);
legend('��ƹ켣')
legend('boxoff')
xlabel('��̷���/km');
ylabel('��չ����/km');
zlabel('�߶ȷ���/km');
title('��ά�ռ��еı�ƹ켣');
grid on