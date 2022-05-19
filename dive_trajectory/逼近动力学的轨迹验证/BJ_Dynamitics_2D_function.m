%-------------------------------------------------------------------------%
%--------------------------��ά����ѧģ��΢�ַ���--------------------------%
% 2022/4/13 16:53 ʷʫ��bug��Lf_nor������L_nor/v1
% Ҫô�ĳ�L_nor/v1��Ҫô����qf
%-------------------------------------------------------------------------%
function state1_dot = BJ_Dynamitics_2D_function(t,state1)
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor
global t0 t0_disp %Ϊ�˹۲�������е�ʲô�ز�

y = state1(1);
h = state1(2);
v = state1(3);
gama = state1(4);

% get_u��ʹ�÷ǽ���ֹ������������ֵ���ڲ�ѯ������ֵ���ڸ�ά���ڵ�����㴦��ֵ�����β�ֵ��
% �Ѻ������ϡ�裬��������߲�ֵ���ȣ�û�仯
t_org = t_nor * (sqrt(R0/g));  %��չ���ʱ������,������һ��
t0 = t * (sqrt(R0/g));         %��չ��ǰʱ��
alpha= interp1(t_org,alpha_nor,t0,'spline');

%���н��ȿ��ӻ�
if (t0 < 1)
    t0_disp = 1;
else
    if (t0 > 1 && t0 > t0_disp)
        fprintf('t0 calculating %f second \n',t0_disp);
        t0_disp = t0_disp + 1;
    end
end

%%
%dynamics f������ѧ
S   = 0.5026;                       %�ο����
% rou = 1.225 * exp(-x(2)/7110);    %�ܶ�rou����Ҫ����ʽ�ƽ�
hf  = -h*R0/7110;                    %h fake���м����
rou = 1.225 * (0.996838143712596 + 0.957192272404239*hf + 0.403293676867969*hf^2 + 0.083714145730035*hf^3 + 0.006865365386321*hf^4);
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
M   = v*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                   %����

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %����

%handle function����������������б�ѩ�����ʽ���бƽ�
sin_x = @(xx) 0.985400513847416*xx - 0.142496853019355*xx^3;
cos_x = @(xx) 0.999396553656190 - 0.495559134405119*xx^2 + 0.036782872656059*xx^4;

sin_g = sin_x(gama);        
cos_g = cos_x(gama);         %gama

%�б�ѩ��ƽ�1/V���Թ�һ��֮���V������[3,7]
V_division = 1.101028089323823 - 0.473922941895471*v + 0.099720451270527*v^2 - 0.010265506659430*v^3 + 4.140771839924636e-04*v^4;

f1 = -v * cos_g;
f2 = v * sin_g;
f3 = -D_nor - sin_g;
f4 = Lf_nor - cos_g*V_division;

state1_dot = zeros(4,1);

%�������t_nor�ڷ���ģ�ΪʲôҪ�ڳ����в��ֵػָ���ʵʱ�䣬�Ѿ�����һ��������
% state1_dot(1) = f1*(sqrt(R0/g));
% state1_dot(2) = f2*(sqrt(R0/g));
% state1_dot(3) = f3*(sqrt(R0/g));
% state1_dot(4) = f4*(sqrt(R0/g));

state1_dot(1) = f1;
state1_dot(2) = f2;
state1_dot(3) = f3;
state1_dot(4) = f4;
end
