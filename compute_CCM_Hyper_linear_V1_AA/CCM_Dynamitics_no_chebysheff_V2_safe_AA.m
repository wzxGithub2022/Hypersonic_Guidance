%% solve local linearization Dynamitics function，解决局部线性化动力学的问题
% 求导会导致逼近效果变差，暂时标记，no chebysheff可以相当改善
% ddf得搞三维矩阵，将雅克比进行到底
% 离谱的alpha bug已修改
% 动力学按照保守写法

clear;clc;
%%
global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

syms h v gama alpha           %暂定为归一化的量，写归一化的动力学

%dynamics f，动力学
S   = 0.5026;                       %参考面积
rou = 1.225 * exp(-h*R0/7110);      %密度rou，可以不用多项式逼近
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %动压
qf  = 0.5 * rou * v*(R0*g);               %q fake，伪动压，已约减速度v，bug标注
M   = v*sqrt(R0*g) / 340;                     %马赫数
m   = 1600;                                   %质量

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor = q*CL*S / (m*g);                        %升力
Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %阻力

%2022/4/6 9:27 OMG，不是吧阿Sir，不会是这里写的bug吧，离大谱
f1 = -v * cos(gama);
f2 = v * sin(gama);
f3 = -D_nor - sin(gama);
f4 = Lf_nor - cos(gama)/v;
f5 = 0;

f_mat(h, v, gama, alpha) = [f1;f2;f3;f4;f5];               %动力学，5*1

b1 = diff(f1,alpha);
b2 = diff(f2,alpha);
b3 = diff(f3,alpha);
b4 = diff(f4,alpha);

B_mat(h, v, gama, alpha) = [b1;b2;b3;b4];

%% df家族
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
df3_alpha = diff(f3,alpha);     %变成b3了

df3 = [df3_y, df3_h, df3_v, df3_gama, df3_alpha];

df4_y = 0;
df4_h = diff(f4,h);
df4_v = diff(f4,v);
df4_gama = diff(f4,gama);
df4_alpha = diff(f4,alpha);     %变成b4了

df4 = [df4_y, df4_h, df4_v, df4_gama, df4_alpha];

df5 = [0,0,0,0,0];

df_mat(h, v, gama, alpha) = [df1;df2;df3;df4;df5];

%% Current state load，导入状态点的信息

load('Trajectory_normalization.mat');

df_mat_value = zeros(5,5,46);
B_mat_value = zeros(4,1,46);    %对AA来说没必要算

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
% save('CCM_Dynamitics_B.mat','B_mat_value');       %对AA来说没必要算
