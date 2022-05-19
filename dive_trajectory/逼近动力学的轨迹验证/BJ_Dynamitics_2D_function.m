%-------------------------------------------------------------------------%
%--------------------------二维动力学模型微分方程--------------------------%
% 2022/4/13 16:53 史诗级bug：Lf_nor不等于L_nor/v1
% 要么改成L_nor/v1，要么修正qf
%-------------------------------------------------------------------------%
function state1_dot = BJ_Dynamitics_2D_function(t,state1)
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor
global t0 t0_disp %为了观察程序运行到什么地步

y = state1(1);
h = state1(2);
v = state1(3);
gama = state1(4);

% get_u，使用非结终止条件的样条插值。在查询点插入的值基于各维中邻点网格点处数值的三次插值。
% 把横坐标变稀疏，或许能提高插值精度，没变化
t_org = t_nor * (sqrt(R0/g));  %扩展标称时间序列,即反归一化
t0 = t * (sqrt(R0/g));         %扩展当前时间
alpha= interp1(t_org,alpha_nor,t0,'spline');

%运行进度可视化
if (t0 < 1)
    t0_disp = 1;
else
    if (t0 > 1 && t0 > t0_disp)
        fprintf('t0 calculating %f second \n',t0_disp);
        t0_disp = t0_disp + 1;
    end
end

%%
%dynamics f，动力学
S   = 0.5026;                       %参考面积
% rou = 1.225 * exp(-x(2)/7110);    %密度rou，需要多项式逼近
hf  = -h*R0/7110;                    %h fake，中间变量
rou = 1.225 * (0.996838143712596 + 0.957192272404239*hf + 0.403293676867969*hf^2 + 0.083714145730035*hf^3 + 0.006865365386321*hf^4);
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %动压
qf  = 0.5 * rou * v*(R0*g);               %q fake，伪动压，已约减速度v，bug标注
M   = v*sqrt(R0*g) / 340;                     %马赫数
m   = 1600;                                   %质量

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor = q*CL*S / (m*g);                        %升力
Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %阻力

%handle function，句柄函数，生成切比雪夫多项式进行逼近
sin_x = @(xx) 0.985400513847416*xx - 0.142496853019355*xx^3;
cos_x = @(xx) 0.999396553656190 - 0.495559134405119*xx^2 + 0.036782872656059*xx^4;

sin_g = sin_x(gama);        
cos_g = cos_x(gama);         %gama

%切比雪夫逼近1/V，对归一化之后的V，区间[3,7]
V_division = 1.101028089323823 - 0.473922941895471*v + 0.099720451270527*v^2 - 0.010265506659430*v^3 + 4.140771839924636e-04*v^4;

f1 = -v * cos_g;
f2 = v * sin_g;
f3 = -D_nor - sin_g;
f4 = Lf_nor - cos_g*V_division;

state1_dot = zeros(4,1);

%你就是拿t_nor在仿真的，为什么要在程序中部分地恢复真实时间，已经换了一个世界了
% state1_dot(1) = f1*(sqrt(R0/g));
% state1_dot(2) = f2*(sqrt(R0/g));
% state1_dot(3) = f3*(sqrt(R0/g));
% state1_dot(4) = f4*(sqrt(R0/g));

state1_dot(1) = f1;
state1_dot(2) = f2;
state1_dot(3) = f3;
state1_dot(4) = f4;
end
