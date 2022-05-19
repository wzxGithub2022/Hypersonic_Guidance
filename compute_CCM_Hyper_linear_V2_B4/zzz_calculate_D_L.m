% 验证归一化动力学的计算结果

global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

load('Trajectory_normalization.mat')
% syms h v gama alpha           %暂定为归一化的量，写归一化的动力学
L_nor = zeros(46,1);
D_nor = zeros(46,1);

for i = 1:46
    
h = Trajectory_normalization(i,3);
v = Trajectory_normalization(i,4);
gama = Trajectory_normalization(i,5);
alpha = Trajectory_normalization(i,6);

%dynamics f，动力学
S   = 0.5026;                       %参考面积
rou = 1.225 * exp(-h*R0/7110);      %密度rou，可以不用多项式逼近
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %动压
qf  = 0.5 * rou * v*(R0*g);               %q fake，伪动压，已约减速度v，bug标注
M   = v*sqrt(R0*g) / 340;                     %马赫数
m   = 1600;                                   %质量

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor(i) = q*CL*S / (m*g);                        %升力
Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor(i) = q*Cd*S / (m*g);                        %阻力

end

t = Trajectory_normalization(:,1) * sqrt(R0/g);

figure(1)
subplot(121),plot(t,L_nor);
title('归一化升力');
subplot(122),plot(t,D_nor);
title('归一化阻力');