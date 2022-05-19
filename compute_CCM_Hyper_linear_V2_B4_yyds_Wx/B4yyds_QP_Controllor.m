%为了提高程序的可比性和运行效率，基于在线算法框架，对B4_fix进行重构
function [k] = B4yyds_QP_Controllor(x_d, u_d, x, lamda, W, fun_f, fun_b)

m = length(x_d);
n = length(u_d);

x_e = x - x_d;

if norm(x_e) < 1e-12         %不能拿一个矢量和一个数字比大小啊喂
    k = 0;
else
    delta_gamma = x_e;
%     delta_gamma = x_e / norm(x_e);  %0426懂了，看论文，推一推，x_e才是正解
    
    Riemann_energy = x_e' * ( W\eye(m) ) * x_e;
    
    A = 2*delta_gamma'*(W\eye(m))*fun_b;
    b = -2*lamda * Riemann_energy + 2*delta_gamma'*(W\eye(m))*(fun_f*zeros(m,1) + fun_b*u_d) - ...
        2*delta_gamma'*(W\eye(m))*(fun_f*(x-x_d) + fun_b*u_d);
    
    H = eye(n);
    f = zeros(n,1);
    
    k = quadprog(H,f,A,b);
end

end
%% 0.01气动
% 射程偏差km
%     0.1559
% 高度偏差km
%    -0.1679
% 落速m/s
%    1.0745e+03
% 落角/度
%   -65.9308
%% norm(x_e)
% 射程偏差km
%     0.0233
% 高度偏差km
%    -0.0747
% 落速m/s
%    1.0786e+03
% 落角/度
%   -65.3745
%% x_e
% 射程偏差km
%   -8.2047e-04
% 高度偏差km
%    -0.0883
% 落速m/s
%    1.0782e+03
% 落角/度
%   -65.3188