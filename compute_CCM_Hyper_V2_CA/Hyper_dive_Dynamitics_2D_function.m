%-------------------------------------------------------------------------%
%--------------------------二维动力学模型微分方程--------------------------%
% 2022/4/13 16:53 史诗级bug：Lf_nor不等于L_nor/v1
% 要么改成L_nor/v1，要么修正qf
%-------------------------------------------------------------------------%
function state1_dot = Hyper_dive_Dynamitics_2D_function(t,state1)
%%
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor
global t0 t0_disp %为了观察程序运行到什么地步

y1 = state1(1);
h1 = state1(2);
v1 = state1(3);
gama1 = state1(4);
alpha1 = state1(5);

% get_u，使用非结终止条件的样条插值。在查询点插入的值基于各维中邻点网格点处数值的三次插值。
% 把横坐标变稀疏，或许能提高插值精度，没变化
t_org = t_nor * (sqrt(R0/g));  %扩展标称时间序列,即反归一化
t0 = t * (sqrt(R0/g));         %扩展当前时间
alpha_std= interp1(t_org,alpha_nor,t0,'spline');

%运行进度可视化
if (t0 < 1)
    t0_disp = 1;
else
    if (t0 > 1 && t0 > t0_disp)
        fprintf('t0 calculating %f second \n',t0_disp);
        t0_disp = t0_disp + 1;
    end
end

% 标称控制留在这个动力学里，反馈控制写在QP里
alpha_QP = alpha_QP_function(t,state1,alpha_std);
% alpha_QP = 0;
alpha_pi = alpha_std + alpha_QP;
%%
%dynamics f，加控制后的动力学
S   = 0.5026;                       %参考面积
rou = 1.225 * exp(-h1*R0/7110);      %密度rou
q   = 0.5 * rou * (v1*sqrt(R0*g))^2;           %动压
qf  = 0.5 * rou * v1*(R0*g);               %q fake，伪动压，已约减速度v，史诗级bug，约减v和约减v1的问题，修正qf
M   = v1*sqrt(R0*g) / 340;                     %马赫数
m   = 1600;                                   %质量

CL  = 0.4172 + 19.41*alpha_pi + 10.17*alpha_pi^2 - M*(0.1004 + 0.7536*alpha_pi);
L_nor = q*CL*S / (m*g);                        %升力
Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %阻力

%% 扰动，先按照你生成CCW的这一套来写，Bw相对独立，易更改
% D1_dtb = D_nor * 0.01 * randn(1)/3;       % randn是正态分布，更真实，但不是非常符合我们的问题
% Lf1_dtb = Lf_nor * 0.01 * randn(1)/3;     %为什么偏差加不上去的样子,细密的步长把随机的效果抵消完了
% D1_dtb = D_nor * 0.01 * (rand(1)-0.5);
% Lf1_dtb = Lf_nor * 0.01 * (rand(1)-0.5);  %尝试均布的干扰，好家伙，看起来有正有负就能消
D1_dtb = D_nor * 0.01;                    %常值干扰，0.01气动力也就(157.8m,167.9m)左右，加控制看看
Lf1_dtb = Lf_nor * 0.01;
% D1_dtb = 0;                               %干扰置零
% Lf1_dtb = 0;    
%%
state1_dot=zeros(5,1);
state1_dot(1) = -v1 * cos(gama1);
state1_dot(2) = v1 * sin(gama1);
state1_dot(3) = -D_nor - sin(gama1) + D1_dtb ;
state1_dot(4) = Lf_nor - cos(gama1)/v1 + Lf1_dtb ;
global step
state1_dot(5) = (alpha_pi - alpha1)/step; %这个写法还是比较科学的，能做到让状态的alpha1

end

%% 看看这样仿真的误差如何，这误差挺大，不是很行，这个值得单独验证和探究一下，BJ逼近动力学的问题已经处理
% % 警告: 在 t=7.368434e-02 处失败。在时间 t 处，若不将步长降至允许的最小值(2.220446e-16)以下，积分公差要求无法满足。 
% > In ode45 (line 360)
%   In Hyper_dive_Trajectory_Simulation (line 29) 

% %dynamics f，动力学
% S   = 0.5026;                       %参考面积
% % rou = 1.225 * exp(-x(2)/7110);    %密度rou，需要多项式逼近
% hf  = h*R0/7110;                    %h fake，中间变量
% rou = 1.225 * (1.124508748184077 + 0.459262515874491*hf + 3.155718757366057*hf^2 + 1.426742501850045*hf^3 + 0.374836249394692*hf^4);
% q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %动压
% qf  = 0.5 * rou * v*(R0*g);               %q fake，伪动压，已约减速度v，bug标注
% M   = v*sqrt(R0*g) / 340;                     %马赫数
% m   = 1600;                                   %质量
% 
% CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
% L_nor = q*CL*S / (m*g);                        %升力
% Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v
% 
% Cd  = 0.3042 + 0.02988*CL^2;
% D_nor = q*Cd*S / (m*g);                        %阻力
% 
% %handle function，句柄函数，生成切比雪夫多项式进行逼近
% sin_x = @(xx) 0.985400513847416*xx - 0.142496853019355*xx^3;
% cos_x = @(xx) 0.999396553656190 - 0.495559134405119*xx^2 + 0.036782872656059*xx^4;
% 
% sin_g = sin_x(gama);        
% cos_g = cos_x(gama);         %gama
% 
% %切比雪夫逼近1/V，对归一化之后的V，区间[3,7]
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