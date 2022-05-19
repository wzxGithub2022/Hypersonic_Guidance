%基于在线算法框架，对CA进行重构
function [k] = CAyyds_QP_Controllor(X_d, u_d, X, alphad, lambda, W)
%% 获得当前状态和标称状态
global R0 g
%alpha输入就是alphad，把标称当成当前攻角，这没问题
%真实状态
y1 = X(1);
h1 = X(2);
v1 = X(3);
gama1 = X(4);
alpha1 = alphad;
%标称状态
y_std = X_d(1);
h_std = X_d(2);
v_std = X_d(3);
gama_std = X_d(4);
alpha_std = X_d(5);
%W_to_M
m = length(X);
n = length(u_d);
CCM_std = W\eye(m);
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
%% 现在可以考虑写成QP的形式了
fx_minus = double( f_mat(h1, v1, gama1, alpha1) ) - double( f_mat(h_std, v_std, gama_std, alpha_std) );
B = double( B_mat(h1, v1, gama1, alpha1) );

x_e = X - X_d(1:4);

if norm(x_e) < 1e-12         %不能拿一个矢量和一个数字比大小啊喂
    k = 0;
else
    delta_gama = x_e;
    Riemann_energy = x_e' * CCM_std * x_e;
    %最小化目标函数
    H = eye(n);
    f = zeros(n,1);
    %不等式约束
    A = 2*delta_gama' * CCM_std * B;
    b = -2*lambda * Riemann_energy - 2*delta_gama' * CCM_std * fx_minus;
    %没有等式约束
    Aeq = [];
    beq = [];
    %没有初值猜测
    lb = [];
    ub = [];
    
    k = quadprog(H,f,A,b,Aeq,beq,lb,ub);
end


end