function [k] = QP_Tracking_Controller(x_e,f_x,M,B,lambda)

delta_gama = x_e / norm(x_e);

%��С��Ŀ�꺯��
H = 2;
f = 0;

%����ʽԼ��
A = delta_gama' * M * B;
b = -lambda * x_e' * M * x_e - delta_gama' * M * f_x * x_e;

%û�е�ʽԼ��
Aeq = [];
beq = [];
%û�г�ֵ�²�
lb = [];
ub = [];

k = quadprog(H,f,A,b,Aeq,beq,lb,ub);

end