%���������㷨��ܣ���CA�����ع�
function [k] = CAyyds_QP_Controllor(X_d, u_d, X, alphad, lambda, W)
%% ��õ�ǰ״̬�ͱ��״̬
global R0 g
%alpha�������alphad���ѱ�Ƶ��ɵ�ǰ���ǣ���û����
%��ʵ״̬
y1 = X(1);
h1 = X(2);
v1 = X(3);
gama1 = X(4);
alpha1 = alphad;
%���״̬
y_std = X_d(1);
h_std = X_d(2);
v_std = X_d(3);
gama_std = X_d(4);
alpha_std = X_d(5);
%W_to_M
m = length(X);
n = length(u_d);
CCM_std = W\eye(m);
%% �Ȳ�����ʡ��һ�������ٶ���
syms h v gama alphaa           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ

%dynamics f������ѧ
S   = 0.5026;                       %�ο����
rou = 1.225 * exp(-h*R0/7110);      %�ܶ�rou�����Բ��ö���ʽ�ƽ�
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
M   = v*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                   %����

CL  = 0.4172 + 19.41*alphaa + 10.17*alphaa^2 - M*(0.1004 + 0.7536*alphaa);
L_nor = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %����

f1 = -v * cos(gama);
f2 = v * sin(gama);
f3 = -D_nor - sin(gama);
f4 = Lf_nor - cos(gama)/v;

f_mat(h, v, gama, alphaa) = [f1;f2;f3;f4];           %����ѧ��4*1����Ҫ�����ã������

b1 = diff(f1,alphaa);
b2 = diff(f2,alphaa);
b3 = diff(f3,alphaa);
b4 = diff(f4,alphaa);

B_mat(h, v, gama, alphaa) = [b1;b2;b3;b4];           % x_dot = Ax+Bu ��B
%% ���ڿ��Կ���д��QP����ʽ��
fx_minus = double( f_mat(h1, v1, gama1, alpha1) ) - double( f_mat(h_std, v_std, gama_std, alpha_std) );
B = double( B_mat(h1, v1, gama1, alpha1) );

x_e = X - X_d(1:4);

if norm(x_e) < 1e-12         %������һ��ʸ����һ�����ֱȴ�С��ι
    k = 0;
else
    delta_gama = x_e;
    Riemann_energy = x_e' * CCM_std * x_e;
    %��С��Ŀ�꺯��
    H = eye(n);
    f = zeros(n,1);
    %����ʽԼ��
    A = 2*delta_gama' * CCM_std * B;
    b = -2*lambda * Riemann_energy - 2*delta_gama' * CCM_std * fx_minus;
    %û�е�ʽԼ��
    Aeq = [];
    beq = [];
    %û�г�ֵ�²�
    lb = [];
    ub = [];
    
    k = quadprog(H,f,A,b,Aeq,beq,lb,ub);
end


end