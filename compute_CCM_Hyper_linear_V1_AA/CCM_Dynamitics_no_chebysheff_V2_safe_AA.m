%% solve local linearization Dynamitics function������ֲ����Ի�����ѧ������
% �󵼻ᵼ�±ƽ�Ч������ʱ��ǣ�no chebysheff�����൱����
% ddf�ø���ά���󣬽��ſ˱Ƚ��е���
% ���׵�alpha bug���޸�
% ����ѧ���ձ���д��

clear;clc;
%%
global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

syms h v gama alpha           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ

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

%2022/4/6 9:27 OMG�����ǰɰ�Sir������������д��bug�ɣ������
f1 = -v * cos(gama);
f2 = v * sin(gama);
f3 = -D_nor - sin(gama);
f4 = Lf_nor - cos(gama)/v;
f5 = 0;

f_mat(h, v, gama, alpha) = [f1;f2;f3;f4;f5];               %����ѧ��5*1

b1 = diff(f1,alpha);
b2 = diff(f2,alpha);
b3 = diff(f3,alpha);
b4 = diff(f4,alpha);

B_mat(h, v, gama, alpha) = [b1;b2;b3;b4];

%% df����
df1_y = 0;
df1_h = 0;
df1_v = diff(f1,v);
df1_gama = diff(f1,gama);
df1_alpha = 0;

df1 = [df1_y, df1_h, df1_v, df1_gama, df1_alpha];

df2_y = 0;
df2_h = 0;
df2_v = diff(f2,v);
df2_gama = diff(f2,gama);
df2_alpha = 0;

df2 = [df2_y, df2_h, df2_v, df2_gama, df2_alpha];

df3_y = 0;
df3_h = diff(f3,h);
df3_v = diff(f3,v);
df3_gama = diff(f3,gama);
df3_alpha = diff(f3,alpha);     %���b3��

df3 = [df3_y, df3_h, df3_v, df3_gama, df3_alpha];

df4_y = 0;
df4_h = diff(f4,h);
df4_v = diff(f4,v);
df4_gama = diff(f4,gama);
df4_alpha = diff(f4,alpha);     %���b4��

df4 = [df4_y, df4_h, df4_v, df4_gama, df4_alpha];

df5 = [0,0,0,0,0];

df_mat(h, v, gama, alpha) = [df1;df2;df3;df4;df5];

%% Current state load������״̬�����Ϣ

load('Trajectory_normalization.mat');

df_mat_value = zeros(5,5,46);
B_mat_value = zeros(4,1,46);    %��AA��˵û��Ҫ��

for i = 1:46
    
    h1 = Trajectory_normalization(i,3);
    v1 = Trajectory_normalization(i,4);
    gama1 = Trajectory_normalization(i,5);
    alpha1 = Trajectory_normalization(i,6);
    state_base = [h1,v1,gama1,alpha1];
    
    df_mat_value(:,:,i) = double( df_mat(h1,v1,gama1,alpha1) );
    B_mat_value(:,:,i) = double( B_mat(h1,v1,gama1,alpha1) );
    
end

%% save and ouput
save('CCM_Dynamitics_df.mat','df_mat_value');
% save('CCM_Dynamitics_B.mat','B_mat_value');       %��AA��˵û��Ҫ��
