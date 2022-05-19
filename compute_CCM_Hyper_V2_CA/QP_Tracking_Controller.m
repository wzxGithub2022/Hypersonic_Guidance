function [k] = QP_Tracking_Controller(x_e,f_x,M,B,lambda)

delta_gama = x_e / norm(x_e);

%最小化目标函数
H = 2;
f = 0;

%不等式约束
A = delta_gama' * M * B;
b = -lambda * x_e' * M * x_e - delta_gama' * M * f_x * x_e;

%没有等式约束
Aeq = [];
beq = [];
%没有初值猜测
lb = [];
ub = [];

k = quadprog(H,f,A,b,Aeq,beq,lb,ub);

end