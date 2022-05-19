%%
clear;clc;
%%
global R0 g
R0 = 10*10^3;
g = 9.81;
load('Trajectory_normalization.mat')
t_nor = Trajectory_normalization(:,1);                %����ʱ������
y_nor = Trajectory_normalization(:,2);
h_nor = Trajectory_normalization(:,3);
v_nor = Trajectory_normalization(:,4);
gama_nor = Trajectory_normalization(:,5);
alpha_nor = Trajectory_normalization(:,6);                %������u

for i=1:1
    
    h = h_nor(i);
    v = v_nor(i);
    gama = gama_nor(i);
    alpha = alpha_nor(i);
    
    %dynamics f������ѧ
    S   = 0.5026;                       %�ο����
    % rou = 1.225 * exp(-x(2)/7110);    %�ܶ�rou����Ҫ����ʽ�ƽ�
    hf  = -h*R0/7110;                    %h fake���м����
    rou = 1.225 * (0.996838143712596 + 0.957192272404239*hf + 0.403293676867969*hf^2 + 0.083714145730035*hf^3 + 0.006865365386321*hf^4);
    q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
    qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
    M   = v*sqrt(R0*g) / 340;                     %�����
    m   = 1600;                                   %����
    
    CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
    L_nor = q*CL*S / (m*g);                        %����
    Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v
    
    Cd  = 0.3042 + 0.02988*CL^2;
    D_nor = q*Cd*S / (m*g);                        %����
    
    %handle function����������������б�ѩ�����ʽ���бƽ�
    sin_x = @(xx) 0.985400513847416*xx - 0.142496853019355*xx^3;
    cos_x = @(xx) 0.999396553656190 - 0.495559134405119*xx^2 + 0.036782872656059*xx^4;
    
    sin_g = sin_x(gama);
    cos_g = cos_x(gama);         %gama
    
    %�б�ѩ��ƽ�1/V���Թ�һ��֮���V������[3,7]
    V_division = 1.101028089323823 - 0.473922941895471*v + 0.099720451270527*v^2 - 0.010265506659430*v^3 + 4.140771839924636e-04*v^4;
    
    f1_BJ(i) = -v * cos_g;
    f2_BJ(i) = v * sin_g;
    f3_BJ(i) = -D_nor - sin_g;
    f4_BJ(i) = Lf_nor - cos_g*V_division;
    
end