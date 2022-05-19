%% solve local linearization Dynamitics function������ֲ����Ի�����ѧ������
% �󵼻ᵼ�±ƽ�Ч������ʱ���
% ddf�ø���ά���󣬽��ſ˱Ƚ��е���

clear;clc;
%%
global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

%handle function��������������������б�ѩ�����ʽ
sin_x = @(x) 0.985400513847416*x - 0.142496853019355*x^3;
cos_x = @(x) 0.999396553656190 - 0.495559134405119*x^2 + 0.036782872656059*x^4;

syms h v gama alpha           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ

%dynamics f������ѧ
S   = 0.5026;                       %�ο����
% rou = 1.225 * exp(-x(2)/7110);    %�ܶ�rou����Ҫ����ʽ�ƽ�
hf  = h*R0/7110;                    %h fake���м����
rou = 1.225 * (1.124508748184077 + 0.459262515874491*hf + 3.155718757366057*hf^2 + 1.426742501850045*hf^3 + 0.374836249394692*hf^4);
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * v*sqrt(R0*g);               %q fake��α��ѹ����Լ���ٶ�v
M   = v*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                   %����

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %����

%�б�ѩ��ƽ�1/V���Թ�һ��֮���V������[3,7]
V_division = 1.101028089323823 - 0.473922941895471*v + 0.099720451270527*v^2 - 0.010265506659430*v^3 + 4.140771839924636e-04*v^4;

sin_g = sin_x(gama);        
cos_g = cos_x(gama);         %gama
sin_a = sin_x(alpha);
cos_a = cos_x(alpha);        %alpha

%2022/4/6 9:27 OMG�����ǰɰ�Sir������������д��bug�ɣ������
f1 = -v * cos_a;
f2 = v * sin_a;
f3 = -D_nor - sin_g;
f4 = Lf_nor - cos_g*V_division;
f5 = 0;

f_mat(h, v, gama, alpha) = [f1;f2;f3;f4;f5];               %����ѧ��5*1

%% df����
df1_y = 0;
df1_h = 0;
df1_v = -cos_a;
df1_gama = 0;
df1_alpha = v * sin_a;

df1 = [df1_y, df1_h, df1_v, df1_gama, df1_alpha];

df2_y = 0;
df2_h = 0;
df2_v = sin_a;
df2_gama = 0;
df2_alpha = v * cos_a;

df2 = [df2_y, df2_h, df2_v, df2_gama, df2_alpha];

df3_y = 0;
df3_h = diff(f3,h);
df3_v = diff(f3,v);
df3_gama = -cos_g;
df3_alpha = diff(f3,alpha);

df3 = [df3_y, df3_h, df3_v, df3_gama, df3_alpha];

df4_y = 0;
df4_h = diff(f4,h);
df4_v = diff(f4,v);
df4_gama = sin_g*V_division;
df4_alpha = diff(f4,alpha);

df4 = [df4_y, df4_h, df4_v, df4_gama, df4_alpha];

df5 = [0,0,0,0,0];

df_mat(h, v, gama, alpha) = [df1;df2;df3;df4;df5];

%% ddf����
ddf1_y = [0,0,0,0,0];
ddf1_h = [0,0,0,0,0];
ddf1_v = [0,0,0,0,sin_a];
ddf1_gama = [0,0,0,0,0];
ddf1_alpha = [0,0,sin_a,0,v*cos_a];

ddf1_mat(h, v, gama, alpha) = [ddf1_y; ddf1_h; ddf1_v; ddf1_gama; ddf1_alpha];

ddf2_y = [0,0,0,0,0];
ddf2_h = [0,0,0,0,0];
ddf2_v = [0,0,0,0,cos_a];
ddf2_gama = [0,0,0,0,0];
ddf2_alpha = [0,0,cos_a,0,-v*sin_a];

ddf2_mat(h, v, gama, alpha) = [ddf2_y; ddf2_h; ddf2_v; ddf2_gama; ddf2_alpha];

ddf3_y = [0,0,0,0,0];
ddf3_h = [0,diff(df3_h,h),diff(df3_h,v),diff(df3_h,gama),diff(df3_h,alpha)];
ddf3_v = [0,diff(df3_v,h),diff(df3_v,v),diff(df3_v,gama),diff(df3_v,alpha)];
ddf3_gama = [0,diff(df3_gama,h),diff(df3_gama,v),diff(df3_gama,gama),diff(df3_gama,alpha)];
ddf3_alpha = [0,diff(df3_alpha,h),diff(df3_alpha,v),diff(df3_alpha,gama),diff(df3_alpha,alpha)];

ddf3_mat(h, v, gama, alpha) = [ddf3_y; ddf3_h; ddf3_v; ddf3_gama; ddf3_alpha];

ddf4_y = [0,0,0,0,0];
ddf4_h = [0,diff(df4_h,h),diff(df4_h,v),diff(df4_h,gama),diff(df4_h,alpha)];
ddf4_v = [0,diff(df4_v,h),diff(df4_v,v),diff(df4_v,gama),diff(df4_v,alpha)];
ddf4_gama = [0,diff(df4_gama,h),diff(df4_gama,v),diff(df4_gama,gama),diff(df4_gama,alpha)];
ddf4_alpha = [0,diff(df4_alpha,h),diff(df4_alpha,v),diff(df4_alpha,gama),diff(df4_alpha,alpha)];

ddf4_mat(h, v, gama, alpha) = [ddf4_y; ddf4_h; ddf4_v; ddf4_gama; ddf4_alpha];

ddf5_mat = zeros(5,5);

% ddf_mat(h, v, gama, alpha) = [ddf1;ddf2;ddf3;ddf4;ddf5];

%% Current state load������״̬�����Ϣ

load('Trajectory_normalization.mat');

f_mat_value = zeros(5,1,46);
df_mat_value = zeros(5,5,46);
ddf_mat_value = zeros(5,5,5,46);

for i = 1:46
    
    h1 = Trajectory_normalization(i,3);
    v1 = Trajectory_normalization(i,4);
    gama1 = Trajectory_normalization(i,5);
    alpha1 = Trajectory_normalization(i,6);
    state_base = [h1,v1,gama1,alpha1];
    
    f_mat_value(:,:,i) = double( f_mat(h1,v1,gama1,alpha1) );
    df_mat_value(:,:,i) = double( df_mat(h1,v1,gama1,alpha1) );
    ddf_mat_value(:,:,1,i) = double( ddf1_mat(h1,v1,gama1,alpha1) );
    ddf_mat_value(:,:,2,i) = double( ddf2_mat(h1,v1,gama1,alpha1) );
    ddf_mat_value(:,:,3,i) = double( ddf3_mat(h1,v1,gama1,alpha1) );
    ddf_mat_value(:,:,4,i) = double( ddf4_mat(h1,v1,gama1,alpha1) );
    
end

%% save and ouput
save('CCM_Dynamitics_f.mat','f_mat_value');
save('CCM_Dynamitics_df.mat','df_mat_value');
save('CCM_Dynamitics_ddf.mat','ddf_mat_value');
