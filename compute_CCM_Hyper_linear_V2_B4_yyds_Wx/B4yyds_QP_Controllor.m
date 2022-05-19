%Ϊ����߳���Ŀɱ��Ժ�����Ч�ʣ����������㷨��ܣ���B4_fix�����ع�
function [k] = B4yyds_QP_Controllor(x_d, u_d, x, lamda, W, fun_f, fun_b)

m = length(x_d);
n = length(u_d);

x_e = x - x_d;

if norm(x_e) < 1e-12         %������һ��ʸ����һ�����ֱȴ�С��ι
    k = 0;
else
    delta_gamma = x_e;
%     delta_gamma = x_e / norm(x_e);  %0426���ˣ������ģ���һ�ƣ�x_e��������
    
    Riemann_energy = x_e' * ( W\eye(m) ) * x_e;
    
    A = 2*delta_gamma'*(W\eye(m))*fun_b;
    b = -2*lamda * Riemann_energy + 2*delta_gamma'*(W\eye(m))*(fun_f*zeros(m,1) + fun_b*u_d) - ...
        2*delta_gamma'*(W\eye(m))*(fun_f*(x-x_d) + fun_b*u_d);
    
    H = eye(n);
    f = zeros(n,1);
    
    k = quadprog(H,f,A,b);
end

end
%% 0.01����
% ���ƫ��km
%     0.1559
% �߶�ƫ��km
%    -0.1679
% ����m/s
%    1.0745e+03
% ���/��
%   -65.9308
%% norm(x_e)
% ���ƫ��km
%     0.0233
% �߶�ƫ��km
%    -0.0747
% ����m/s
%    1.0786e+03
% ���/��
%   -65.3745
%% x_e
% ���ƫ��km
%   -8.2047e-04
% �߶�ƫ��km
%    -0.0883
% ����m/s
%    1.0782e+03
% ���/��
%   -65.3188