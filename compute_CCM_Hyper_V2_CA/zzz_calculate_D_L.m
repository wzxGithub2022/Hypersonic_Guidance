% 原始动力学的计算结果
% 算一算这升力和阻力到底有多大
global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

load('Trajectory_normalization.mat')
L1 = zeros(46,1);
D1 = zeros(46,1);

for i = 1:46
    
yy = Trajectory_normalization(i,2)*R0;
h = Trajectory_normalization(i,3)*R0;
v = Trajectory_normalization(i,4)*sqrt(R0*g);
gama = Trajectory_normalization(i,5);
alpha = Trajectory_normalization(i,6);

y = [yy,h,v,gama];         %初值航程35km 高度20km 速度 航迹角（弹道倾角)
ui = alpha;

S   = 0.5026;                       %参考面积
rou = 1.225 * exp(-y(2)/7110);    %密度rou
q   = 0.5 * rou * y(3)^2;           %动压    
M   = y(3)/340;                 %马赫数
m   = 1600;                     %质量
CL  = 0.4172 + 19.41*ui + 10.17*ui^2 - M*(0.1004+0.7536*ui);     %u为控制量：u=alpha
L   = q*CL*S;                   %升力
Cd0 = 0.3042;
Cd  = Cd0 + 0.02988*CL^2;         %阻力系数
D   = q*Cd*S;                   %阻力

L1(i) = L / (m*g);     %-3.5276
D1(i) = D / (m*g);     %1.1999，升力和阻力都这么大

end

t = Trajectory_normalization(:,1) * sqrt(R0/g);

figure(1)
subplot(121),plot(t,L1);
title('归一化升力');
subplot(122),plot(t,D1);
title('归一化阻力');