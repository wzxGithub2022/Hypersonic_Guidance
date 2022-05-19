% ��֤��һ������ѧ�ļ�����

global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

load('Trajectory_normalization.mat')
% syms h v gama alpha           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ
L_nor = zeros(46,1);
D_nor = zeros(46,1);

for i = 1:46
    
h = Trajectory_normalization(i,3);
v = Trajectory_normalization(i,4);
gama = Trajectory_normalization(i,5);
alpha = Trajectory_normalization(i,6);

%dynamics f������ѧ
S   = 0.5026;                       %�ο����
rou = 1.225 * exp(-h*R0/7110);      %�ܶ�rou�����Բ��ö���ʽ�ƽ�
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
M   = v*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                   %����

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor(i) = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor(i) = q*Cd*S / (m*g);                        %����

end

t = Trajectory_normalization(:,1) * sqrt(R0/g);

figure(1)
subplot(121),plot(t,L_nor);
title('��һ������');
subplot(122),plot(t,D_nor);
title('��һ������');