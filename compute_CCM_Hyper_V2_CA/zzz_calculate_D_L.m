% ԭʼ����ѧ�ļ�����
% ��һ�������������������ж��
global R0 g
R0 = 10*10^3;                    %R0��λ��m
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

y = [yy,h,v,gama];         %��ֵ����35km �߶�20km �ٶ� �����ǣ��������)
ui = alpha;

S   = 0.5026;                       %�ο����
rou = 1.225 * exp(-y(2)/7110);    %�ܶ�rou
q   = 0.5 * rou * y(3)^2;           %��ѹ    
M   = y(3)/340;                 %�����
m   = 1600;                     %����
CL  = 0.4172 + 19.41*ui + 10.17*ui^2 - M*(0.1004+0.7536*ui);     %uΪ��������u=alpha
L   = q*CL*S;                   %����
Cd0 = 0.3042;
Cd  = Cd0 + 0.02988*CL^2;         %����ϵ��
D   = q*Cd*S;                   %����

L1(i) = L / (m*g);     %-3.5276
D1(i) = D / (m*g);     %1.1999����������������ô��

end

t = Trajectory_normalization(:,1) * sqrt(R0/g);

figure(1)
subplot(121),plot(t,L1);
title('��һ������');
subplot(122),plot(t,D1);
title('��һ������');