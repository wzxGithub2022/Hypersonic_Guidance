function [k] = alpha_QP_function(t,state1,alpha_std)
%%
% ˼�룺
% Ҫ��ƫ���׼_std�͵�ǰ
% Ҫ�ж���ѧ���÷��ź�������ʽд���������Է�������
% ��CCW��CCM�����棬�ٲ�ֵ
% ���룺
% ��ʵ����t ��һ��ʱ�� state1 ��һ��״̬ 
% ϵͳ��������ʱ���״̬(��ֵ)�ͱ�ƵĿ���(alpha_std)
% lambda�ǹ̶�0.5�ܵģ�CCM��ֵ�õ���ǰʱ���CCM
%% ��õ�ǰ״̬�ͱ��״̬
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor
%��ʵ״̬
y1 = state1(1);
h1 = state1(2);
v1 = state1(3);
gama1 = state1(4);
alpha1 = state1(5);
%���״̬
t_org = t_nor * (sqrt(R0/g));  %��չ���ʱ�����У�������һ��
t0 = t * (sqrt(R0/g));         %��չ��ǰ��һ��ʱ�䣬�õ���ǰ��ʵʱ��
y_std= interp1(t_org,y_nor,t0,'spline');
h_std= interp1(t_org,h_nor,t0,'spline');
v_std= interp1(t_org,v_nor,t0,'spline');
gama_std= interp1(t_org,gama_nor,t0,'spline');
% alpha_stdֱ����
CCM_std = zeros(4,4);
for ccm_line = 1:4
    for ccm_row = 1:4
        ccm_ij_temp = reshape( CCM_nor(ccm_line,ccm_row,:),46,1 );  %reshape�Ǹ��ö���
        CCM_std(ccm_line,ccm_row) = interp1(t_org,ccm_ij_temp,t0,'spline');
    end
end
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
%% ���ڿ��Կ���д��QP����ʽ�ˣ���ǰ�Ĺ�����ô��ã���д״̬
fx_minus = double( f_mat(h1, v1, gama1, alpha1) ) - double( f_mat(h_std, v_std, gama_std, alpha_std) );
B = double( B_mat(h1, v1, gama1, alpha1) );

x1 = [y1;h1;v1;gama1];
x_std = [y_std;h_std;v_std;gama_std];
x_e = x1 - x_std;

if norm(x_e) == 0
    k = 0;
else
%     delta_gama = x_e / norm(x_e);   %0426���ˣ������ģ���һ�ƣ�x_e��������
    delta_gama = x_e;
    
    % CCM_std�Ѿ�д����
    lambda = 0.7;
    
    %��С��Ŀ�꺯��
    H = 2;
    f = 0;
    
    %����ʽԼ��
    A = delta_gama' * CCM_std * B;
    b = -lambda * x_e' * CCM_std * x_e - delta_gama' * CCM_std * fx_minus;
    
    %û�е�ʽԼ��
    Aeq = [];
    beq = [];
    %û�г�ֵ�²�
    lb = [];
    ub = [];
    
    k = quadprog(H,f,A,b,Aeq,beq,lb,ub);
end

end