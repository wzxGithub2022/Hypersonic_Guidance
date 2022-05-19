%���������㷨��ܣ���CA�����ع�
function [dX ,alpha] = CAyyds_Dynamics(t,X)     %alpha��д���������ʹ��

global t_r u_r y_r h_r v_r gamma_r feedback disturb
global R0 g

fprintf('time now is %f second \n',t*sqrt(R0/g));

n = 4;
lambda = 0.9;
condn_prev = 11;
ccm_eps = 0.05;
return_metric = 1;

yd = interp1(t_r , y_r ,t,'spline');     %Ĭ�Ϸ��������Բ�ֵ
hd = interp1(t_r , h_r ,t,'spline');
vd = interp1(t_r , v_r ,t,'spline');
gammad = interp1(t_r , gamma_r ,t,'spline');
alphad = interp1(t_r ,u_r, t,'spline');

%��ÿ����ֵ����CCM���������������46��CCM��֮���ٶ�CCM��ֵ����ȷʵ�̻������߼���+����Ӧ�õĿ����
X_d = [yd;hd;vd;gammad;alphad];

if feedback
    u_d = 0;     %��Ƶķ�����û�еģ�ֱ�Ӿ����ù��ǳ��򻯵ظı�
%     [condn] = CAyyds_find_condn(n,X_d,condn_prev,lambda,ccm_eps);
    condn = condn_prev;
    [~,~,~,W] = CAyyds_CCM_Opt(n,X_d,condn,lambda,ccm_eps,return_metric);
    [k] = CAyyds_QP_Controllor(X_d, u_d, X, alphad, lambda, W);
    alpha = alphad+k;
else
    alpha = alphad;
end

y = X(1);
h = X(2);
v = X(3);
gama = X(4);

%dynamics f������ѧ
S   = 0.5026;                       %�ο����
rou = 1.225 * exp(-h*R0/7110);      %�ܶ�rou�����Բ��ö���ʽ�ƽ�
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
M   = v*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                   %����

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %����

if disturb
    D1_dtb = D_nor * 0.01;
    Lf1_dtb = Lf_nor * 0.01;
else
    D1_dtb = 0;
    Lf1_dtb = 0;
end

dX = [-v * cos(gama);
      v * sin(gama);
      -D_nor - sin(gama) + D1_dtb;
      Lf_nor - cos(gama)/v + Lf1_dtb];

end