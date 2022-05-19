%% compare state_dot
% 轨迹验证固然是最可信且可视化的方法，我们先低级一点，比较计算结果嘛
clear;clc;
%%
global R0 g
R0 = 10*10^3;
g = 9.81;
load('Trajectory_normalization.mat')
t_nor = Trajectory_normalization(:,1);                %控制时间序列
y_nor = Trajectory_normalization(:,2);
h_nor = Trajectory_normalization(:,3);
v_nor = Trajectory_normalization(:,4);
gama_nor = Trajectory_normalization(:,5);
alpha_nor = Trajectory_normalization(:,6);                %控制量u

for i=1:46
    
    h = h_nor(i);
    v = v_nor(i);
    gama = gama_nor(i);
    alpha = alpha_nor(i);
    
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
    
    f1_BJ(i) = -v * cos_g;
    f2_BJ(i) = v * sin_g;
    f3_BJ(i) = -D_nor - sin_g;
    f4_BJ(i) = Lf_nor - cos_g*V_division;
    
end

for i=1:46
    
    h1 = h_nor(i);
    v1 = v_nor(i);
    gama1 = gama_nor(i);
    alpha_pi = alpha_nor(i);
    
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
    
    f1(i) = -v1 * cos(gama1);
    f2(i) = v1 * sin(gama1);
    f3(i) = -D_nor - sin(gama1);
    f4(i) = Lf_nor - cos(gama1)/v1;
    
end

for i = 1:46
   f1_err(i) = (f1_BJ(i) - f1(i)) /  f1(i);
   f2_err(i) = (f2_BJ(i) - f2(i)) /  f2(i);
   f3_err(i) = (f3_BJ(i) - f3(i)) /  f3(i);
   f4_err(i) = (f4_BJ(i) - f4(i)) /  f4(i);
end

figure(1)
subplot(221),plot(f1_err);title('f1相对误差');
subplot(222),plot(f2_err);title('f2相对误差');
subplot(223),plot(f3_err);title('f3相对误差');
subplot(224),plot(f4_err);title('f4相对误差');

% 这就是大问题啊
% 改了错误，但这逼近效果固然不怎么样，误差能达到4%