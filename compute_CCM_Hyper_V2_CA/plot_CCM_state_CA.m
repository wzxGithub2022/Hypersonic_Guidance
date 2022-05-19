%%
load('state_CCM_CA.mat');
x = linspace(1,46,46);
%%
% alpha_w = max(sigma_ThBw(:));       %(:)表示其中所有元素并用列表示
% d_bar = alpha_w/lambda;             %这d_bar没有乘w上界，属于是J_CCM也不用再除w上界了，一步到位
% disp('d_bar'); 
% disp(d_bar);
figure(1)
plot(x,state_CCM(1,:));
title('d bar');
%%
% disp('Control:'); 
% disp(max(d_bar*delta_u(:)));        %这是控制上界
figure(2)
plot(x,state_CCM(2,:),'-o','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('选取的第i个状态点');
ylabel('控制上界');
title('参数选取导致的计算奇异现象');
grid on
%%
% disp('W:'); 
% disp(min(min(min(eig_W(:,:,:,1)))));      %小中小
% disp(max(max(max(eig_W(:,:,:,2)))));      %大中大
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
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

t = Trajectory_normalization(:,1)*sqrt(R0/g);

figure(6)
plot(t,state_CCM(1,:),'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('时间t/s');
ylabel('Tube \alpha_w/\lambda');
title('基于单位扰动的Tube');
grid on

figure(7)
subplot(221),
plot(t,state_CCM(7,:)*R0,'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('时间t/s');
ylabel('Tube 射程y/m');
title('基于单位扰动的Tube');
grid on

subplot(222),
plot(t,state_CCM(8,:)*R0,'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('时间t/s');
ylabel('Tube 高度h/m');
title('基于单位扰动的Tube');
grid on

subplot(223),
plot(t,state_CCM(9,:)*sqrt(R0*g),'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('时间t/s');
ylabel('Tube 速度v/(m/s)');
title('基于单位扰动的Tube');
grid on

subplot(224),
plot(t,rad2deg( state_CCM(10,:) ),'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('时间t/s');
ylabel('Tube 弹道倾角\gamma/度');
title('基于单位扰动的Tube');
grid on
%%
figure(8)
plot(t,( state_CCM(7,:) - state_CCM(8,:) )*R0,'-','Color',[0 0.447 0.741],'Linewidth',2);
xlabel('时间t/s');
ylabel('射程与高度的Tube偏差/m');
title('基于单位扰动的Tube');
grid on
%%
clc
x = Trajectory_normalization(:,2)*R0/1000;
y = zeros(length(x));
z = Trajectory_normalization(:,3)*R0/1000;
figure(9)
L0 = plot3(x,y,z,'--*','Color',[0 0.447 0.741],'Linewidth',1);

% legend(L0 ,'标称轨迹')
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

xlabel('射程方向/km');
ylabel('拓展方向/km');
zlabel('高度方向/km');
title('基于标称轨迹的三维扩展Tube');
grid on

% legend(L0 ,'标称轨迹','Location','northeast')
% legend('boxoff')
% ah=axes('position',get(gca,'position'),'visible','off');
% legend(ah,L1 ,'Tube','Location','northeast')
% legend('boxoff')

% legend([L0 L1],{'标称轨迹','Tube'})
%%
clear;clc;close all

load('Trajectory_normalization.mat')

global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

x = Trajectory_normalization(:,2)*R0/1000;
y = zeros(length(x));
z = Trajectory_normalization(:,3)*R0/1000;
figure(10)

plot3(x,y,z,'--*','Color',[0 0.447 0.741],'Linewidth',1);
legend('标称轨迹')
legend('boxoff')
xlabel('射程方向/km');
ylabel('拓展方向/km');
zlabel('高度方向/km');
title('三维空间中的标称轨迹');
grid on