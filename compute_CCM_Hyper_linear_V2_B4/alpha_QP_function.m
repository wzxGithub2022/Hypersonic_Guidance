function [k] = alpha_QP_function(t,state1,alpha_std)
%% B4和CA只有只有QP是不一样的，添加控制器的轨迹仿真那都一样
% 思想：
% 要有偏差，标准_std和当前
% 要有动力学，用符号函数的形式写，这样可以反复调用
% 从CCW到CCM先求逆，再插值
% 输入：
% 真实量：t 归一化时间 state1 归一化状态 
% 系统储存量：时序的状态(插值)和标称的控制(alpha_std)
% lambda是固定0.5跑的，CCM插值得到当前时间的CCM
%% 获得当前状态和标称状态
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor B_nor df_nor
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
B_std = zeros(4,1);
for B_num = 1:4
        B_ij_temp = reshape( B_nor(B_num,1,:),46,1 );  %reshape是个好东西
        B_std(B_num,1) = interp1(t_org,B_ij_temp,t0,'spline');
end
df_std = zeros(4,4);
for df_line = 1:4
    for df_row = 1:4
        df_ij_temp = reshape( df_nor(df_line,df_row,:),46,1 );  %reshape是个好东西
        df_std(ccm_line,ccm_row) = interp1(t_org,df_ij_temp,t0,'spline');
    end
end
%% 现在可以考虑写成QP的形式了，当前的攻角怎么获得，改写状态

x1 = [y1;h1;v1;gama1];
x_std = [y_std;h_std;v_std;gama_std];
x_e = x1 - x_std;

if norm(x_e) == 0
    k = 0;              %没有误差直接这样还快啊
else
    delta_gama = x_e / norm(x_e);   %如果没有误差，norm是0，作为除数，程序无法运行
    
    % CCM_std已经写好了
    lambda = 0.5;
    
    %最小化目标函数
    H = 2;
    f = 0;
    
    %不等式约束
    A = delta_gama' * CCM_std * B_std;
    b = -lambda * x_e' * CCM_std * x_e - delta_gama' * CCM_std * df_std * x_e;
    
    %没有等式约束
    Aeq = [];
    beq = [];
    %没有初值猜测
    lb = [];
    ub = [];
    
    k = quadprog(H,f,A,b,Aeq,beq,lb,ub);
end

end