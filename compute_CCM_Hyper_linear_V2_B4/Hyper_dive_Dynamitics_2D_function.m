%-------------------------------------------------------------------------%
%--------------------------��ά����ѧģ��΢�ַ���--------------------------%
% 2022/4/13 16:53 ʷʫ��bug��Lf_nor������L_nor/v1
% Ҫô�ĳ�L_nor/v1��Ҫô����qf
%-------------------------------------------------------------------------%
function state1_dot = Hyper_dive_Dynamitics_2D_function(t,state1)
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor B_nor df_nor
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
alpha_QP = alpha_QP_function(t,state1,alpha_std);     %�������Ƶ�Ч��
% alpha_QP = 0;                                           %����Ӱ���Ч��
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
D1_dtb = D_nor * 0.01;                    %��ֵ���ţ�0.01������Ҳ��(157.8m,167.9m)���ң��ӿ��ƿ���
Lf1_dtb = Lf_nor * 0.01;
% D1_dtb = 0;                               %��������
% Lf1_dtb = 0;       
%%
state1_dot=zeros(5,1);
state1_dot(1) = -v1 * cos(gama1);
state1_dot(2) = v1 * sin(gama1);
state1_dot(3) = -D_nor - sin(gama1) + D1_dtb ;
state1_dot(4) = Lf_nor - cos(gama1)/v1 + Lf1_dtb ;
global step
state1_dot(5) = (alpha_pi - alpha1)/step; %���д�����ǱȽϿ�ѧ�ģ���������״̬��alpha1

end

%% ������������������Σ������ͦ�󣬲��Ǻ���
% % ����: �� t=7.368434e-02 ��ʧ�ܡ���ʱ�� t �������������������������Сֵ(2.220446e-16)���£����ֹ���Ҫ���޷����㡣 
% > In ode45 (line 360)
%   In Hyper_dive_Trajectory_Simulation (line 29) 

% %dynamics f������ѧ
% S   = 0.5026;                       %�ο����
% % rou = 1.225 * exp(-x(2)/7110);    %�ܶ�rou����Ҫ����ʽ�ƽ�
% hf  = h*R0/7110;                    %h fake���м����
% rou = 1.225 * (1.124508748184077 + 0.459262515874491*hf + 3.155718757366057*hf^2 + 1.426742501850045*hf^3 + 0.374836249394692*hf^4);
% q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
% qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
% M   = v*sqrt(R0*g) / 340;                     %�����
% m   = 1600;                                   %����
% 
% CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
% L_nor = q*CL*S / (m*g);                        %����
% Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v
% 
% Cd  = 0.3042 + 0.02988*CL^2;
% D_nor = q*Cd*S / (m*g);                        %����
% 
% %handle function����������������б�ѩ�����ʽ���бƽ�
% sin_x = @(xx) 0.985400513847416*xx - 0.142496853019355*xx^3;
% cos_x = @(xx) 0.999396553656190 - 0.495559134405119*xx^2 + 0.036782872656059*xx^4;
% 
% sin_g = sin_x(gama);        
% cos_g = cos_x(gama);         %gama
% 
% %�б�ѩ��ƽ�1/V���Թ�һ��֮���V������[3,7]
% V_division = 1.101028089323823 - 0.473922941895471*v + 0.099720451270527*v^2 - 0.010265506659430*v^3 + 4.140771839924636e-04*v^4;
% 
% f1 = -v * cos_g;
% f2 = v * sin_g;
% f3 = -D_nor - sin_g;
% f4 = Lf_nor - cos_g*V_division;
% 
% state1_dot(1) = f1;
% state1_dot(2) = f2;
% state1_dot(3) = f3;
% state1_dot(4) = f4;