%-------------------------------------------------------------------------%
%--------------------------��ά����ѧģ��΢�ַ���--------------------------%
% 2022/4/13 16:53 ʷʫ��bug��Lf_nor������L_nor/v1
% Ҫô�ĳ�L_nor/v1��Ҫô����qf
%-------------------------------------------------------------------------%
function state1_dot = Hyper_dive_Dynamitics_2D_function(t,state1)
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor
global t0 t0_disp %Ϊ�˹۲�������е�ʲô�ز�

y1 = state1(1);
h1 = state1(2);
v1 = state1(3);
gama1 = state1(4);
alpha1 = state1(5);

% get_u��ʹ�÷ǽ���ֹ������������ֵ���ڲ�ѯ������ֵ���ڸ�ά���ڵ�����㴦��ֵ�����β�ֵ��
% �Ѻ������ϡ�裬��������߲�ֵ���ȣ�û�仯
t_org = t_nor * (sqrt(R0/g));  %��չ���ʱ������,������һ��
t0 = t * (sqrt(R0/g));         %��չ��ǰʱ��
alpha_std= interp1(t_org,alpha_nor,t0,'spline');

%���н��ȿ��ӻ�
if (t0 < 1)
    t0_disp = 1;
else
    if (t0 > 1 && t0 > t0_disp)
        fprintf('t0 calculating %f second \n',t0_disp);
        t0_disp = t0_disp + 1;
    end
end

% ��ƿ��������������ѧ���������д��QP��
alpha_QP = 0;
alpha_pi = alpha_std + alpha_QP;
%%
%dynamics f���ӿ��ƺ�Ķ���ѧ
S   = 0.5026;                       %�ο����
rou = 1.225 * exp(-h1*R0/7110);      %�ܶ�rou
q   = 0.5 * rou * (v1*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * v1*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��ʷʫ��bug��Լ��v��Լ��v1�����⣬����qf
M   = v1*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                   %����

CL  = 0.4172 + 19.41*alpha_pi + 10.17*alpha_pi^2 - M*(0.1004 + 0.7536*alpha_pi);
L_nor = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %����

%% �Ŷ����Ȱ���������CCW����һ����д��Bw��Զ������׸���
% D1_dtb = D_nor * 0.01 * randn(1)/3;       % randn����̬�ֲ�������ʵ�������Ƿǳ��������ǵ�����
% Lf1_dtb = Lf_nor * 0.01 * randn(1)/3;     %Ϊʲôƫ��Ӳ���ȥ������,ϸ�ܵĲ����������Ч����������
% D1_dtb = D_nor * 0.01 * (rand(1)-0.5);
% Lf1_dtb = Lf_nor * 0.01 * (rand(1)-0.5);  %���Ծ����ĸ��ţ��üһ�����������и�������
% D1_dtb = D_nor * 0.01;                    %��ֵ���ţ�0.01������Ҳ��(157.8m,167.9m)���ң��ӿ��ƿ���
% Lf1_dtb = Lf_nor * 0.01;
D1_dtb = 0;                               %��������
Lf1_dtb = 0;    
%%
state1_dot=zeros(5,1);
state1_dot(1) = -v1 * cos(gama1);
state1_dot(2) = v1 * sin(gama1);
state1_dot(3) = -D_nor - sin(gama1) + D1_dtb ;
state1_dot(4) = Lf_nor - cos(gama1)/v1 + Lf1_dtb ;
global step
state1_dot(5) = (alpha_pi - alpha1)/step; %���д�����ǱȽϿ�ѧ�ģ���������״̬��alpha1

end