function [k] = alpha_QP_function(t,state1,alpha_std)
%%
% 思想：
% 要有偏差，标准_std和当前
% 要有动力学，用符号函数的形式写，这样可以反复调用
% 从CCW到CCM先求逆，再插值
% 输入：
% 真实量：t 归一化时间 state1 归一化状态 
% 系统储存量：时序的状态(插值)和标称的控制(alpha_std)
% lambda是固定0.5跑的，CCM插值得到当前时间的CCM
%% 获得当前状态和标称状态
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor
%真实状态
y1 = state1(1);
h1 = state1(2);
v1 = state1(3);
gama1 = state1(4);
alpha1 = state1(5);
%标称状态
t_org = t_nor * (sqrt(R0/g));  %扩展标称时间序列，即反归一化
t0 = t * (sqrt(R0/g));         %扩展当前归一化时间，得到当前真实时间
y_std= interp1(t_org,y_nor,t0,'spline');
h_std= interp1(t_org,h_nor,t0,'spline');
v_std= interp1(t_org,v_nor,t0,'spline');
gama_std= interp1(t_org,gama_nor,t0,'spline');
% alpha_std直接用
CCM_std = zeros(4,4);
for ccm_line = 1:4
    for ccm_row = 1:4
        ccm_ij_temp = reshape( CCM_nor(ccm_line,ccm_row,:),46,1 );  %reshape是个好东西
        CCM_std(ccm_line,ccm_row) = interp1(t_org,ccm_ij_temp,t0,'spline');
    end
end
%% 先不考虑省这一点运行速度了
syms h v gama alphaa           %暂定为归一化的量，写归一化的动力学

%dynamics f，动力学
S   = 0.5026;                       %参考面积
rou = 1.225 * exp(-h*R0/7110);      %密度rou，可以不用多项式逼近
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %动压
qf  = 0.5 * rou * v*(R0*g);               %q fake，伪动压，已约减速度v，bug标注
M   = v*sqrt(R0*g) / 340;                     %马赫数
m   = 1600;                                   %质量

CL  = 0.4172 + 19.41*alphaa + 10.17*alphaa^2 - M*(0.1004 + 0.7536*alphaa);
L_nor = q*CL*S / (m*g);                        %升力
Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %阻力

f1 = -v * cos(gama);
f2 = v * sin(gama);
f3 = -D_nor - sin(gama);
f4 = Lf_nor - cos(gama)/v;

f_mat(h, v, gama, alphaa) = [f1;f2;f3;f4];           %动力学，4*1，但要他无用，不输出

b1 = diff(f1,alphaa);
b2 = diff(f2,alphaa);
b3 = diff(f3,alphaa);
b4 = diff(f4,alphaa);

B_mat(h, v, gama, alphaa) = [b1;b2;b3;b4];           % x_dot = Ax+Bu 的B
%% 现在可以考虑写成QP的形式了，当前的攻角怎么获得，改写状态
fx_minus = double( f_mat(h1, v1, gama1, alpha1) ) - double( f_mat(h_std, v_std, gama_std, alpha_std) );
B = double( B_mat(h1, v1, gama1, alpha1) );

x1 = [y1;h1;v1;gama1];
x_std = [y_std;h_std;v_std;gama_std];
x_e = x1 - x_std;

if norm(x_e) == 0
    k = 0;
else
%     delta_gama = x_e / norm(x_e);   %0426懂了，看论文，推一推，x_e才是正解
    delta_gama = x_e;
    
    % CCM_std已经写好了
    lambda = 0.7;
    
    %最小化目标函数
    H = 2;
    f = 0;
    
    %不等式约束
    A = delta_gama' * CCM_std * B;
    b = -lambda * x_e' * CCM_std * x_e - delta_gama' * CCM_std * fx_minus;
    
    %没有等式约束
    Aeq = [];
    beq = [];
    %没有初值猜测
    lb = [];
    ub = [];
    
    k = quadprog(H,f,A,b,Aeq,beq,lb,ub);
end

end