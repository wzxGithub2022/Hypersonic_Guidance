%% 不归一化会算力障碍，弃置，留底
% 2022/3/30 19:21

%% compute_CCM_Hyper
% 主程序里完全没有用到spotless的东西，他们只是优化用的工具
clear all; close all; clc;
% yalmip('clear'); %画图包里的功能，暂且不用
warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load，导入状态点的信息

load('Trajectory_all_information.mat');

i = 1;
h0 = Trajectory_all_information(i,3);
v0 = Trajectory_all_information(i,4);
gama0 = Trajectory_all_information(i,5);
alpha0 = Trajectory_all_information(i,6);
state_base = [h0,v0,gama0,alpha0];

%% Constants，常数

n = 5;                      %动力学有5个量
g = 9.81;
ccm_eps = 0.05;

%注意，这是每一点的误差限
h_lim = 1;                  %高度h单位：km
v_lim = 100;
gama_lim = pi/18;           %这两个角度的单位都是：弧度
alpha_lim = pi/36;

%% Uncomment for global optimization of bound，这就是OPT-CCM

%寻找最小条件数的过程，先1.2倍迭代增大，再二分法搜索
%优化的终点是什么，在不同的收缩率下，CCM条件数达标时，性能指标是不同的，给一系列lambda，画出对应的性能指标，人工挑选最优的lambda
%euc前期是sqrt(double(w_upper/w_lower))/lambda，等价于（1/lambda^2）*(w_upper/w_lower)，后期是d_bar * diag(W_upper)^1/2
%d_bars前期是sqrt(double(1/w_lower))/lambda，后期是alpha_w/lambda

lambda_range = linspace(0.7,0.95,5);              %线性插值，lambda的范围
lambda_range = (1/100)*round(lambda_range*100);   %四舍五入为整数
euc_bounds = NaN(length(lambda_range),1);         %注意，前期优化寻找CCM和后期准确计算时的含义不严格相同，
d_bars = NaN(length(lambda_range),1);             %bar是上面一横
cond_bound = NaN(length(lambda_range),1);         %条件数

eps = 1;
condn_prev = 50;                                  %预设条件数
return_metric = 0;

for ll = 1:length(lambda_range)       %对每一个lambda操作，是小写的“LL”,美观，接受了
    
    lambda = lambda_range(ll);
    fprintf('**********\n');
    fprintf('lambda: %f\n', lambda);
    solved = 0;
    
    %Determine upper bound，决定上界
    cond_l = condn_prev;              
    cond_u = 1.2*condn_prev;
%------------------------若代码一直在循环就是卡在这个区间，sos问题一直得不到解决-----------------------%
    while (~solved) 
        fprintf(' cond_u: %.2f: ', cond_u);
        [sos_prob,~,~] = CCM_Hyper_Opt(n,g,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                cond_u,lambda,ccm_eps,return_metric);
        if (sos_prob == 0)
            solved = 1;
            fprintf('feasible \n');
        else
            %shift up condition number range
            %解决不了的话就一点点调大condition number，只调大不一定解决问题
            fprintf('你看这解不出来了吧\n');
            cond_l = cond_u;
            cond_u = 1.2*cond_u;
        end
    end
%---------------------------------------------------------------------------------------------------%
    if (solved)           %直到解决了，把euc_bounds(ll)定下来
        euc_bounds(ll) = sqrt(cond_u)/lambda;       %初步算出来一个值，此时还无法计算d_bars
        fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
    else
        continue;
    end
    
    %Now do bisection search，做二分搜索
    while(cond_u - cond_l >= eps)     %误差较大，细分优化
        condn = (cond_l+cond_u)/2;
        fprintf(' cond: %.2f', condn);
        [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,g,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                condn,lambda,ccm_eps,return_metric);
        if (sos_prob == 0)
            fprintf(' feasible\n');       %这就可行了，向左再取小一点，向左动的话，euc和d就还能再优化一些      
            euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;   %euc等价于（1/lambda^2）*(w_upper/w_lower)
            d_bars(ll) = sqrt(double(1/w_lower))/lambda;             %这里才首次计算d_bars
            cond_u = condn;
        else
            fprintf(' infeasible\n');     %不可行的话，向右放大一点
            cond_l = condn;
        end
    end

    condn_prev = cond_u;              %误差较小，全都相等
    cond_bound(ll) = cond_u;

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end

pause; 

%% Pick a (lambda,condn)，这是组合拳，挑一套的

lambda = 0.83;          %这是挑出来的
condn = 132.8;
return_metric = 1;

save_name = 'metric_Hyper_vectorized.mat';          %之后的部分只调用了一次CCM_OTP，还用save_name把东西取出来，之前打掉的部分已经完成了优化CCM的功能
[sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,g,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                condn,lambda,ccm_eps,return_metric,save_name);      %在这一步存好，sos_prob:sum of squares problem

%% Compute aux control bound，辅助控制边界
%这是想看看在算出来的(lambda，condn)情况下，CCM有什么样的表现，需要多大的u才能控得住
%控制和扰动是耦合的，给你一个扰动上界，你还想把它控制在什么样tube内，决定了你的u的情况

load(save_name);
disp('Checking CCM conditions and Computing control bound...');
lambda = 0.998*lambda;      %缩小一丢丢,这不是很有意义

B = [zeros(4,1);1];       % 5*1 的double
 
Bw = @(x)[zeros(2,2);
           -1/1600,0;
     0,1/(1600*2000);
                0,0];      %这样的写法直接就是小function啊，气动力扰动对系统的输入，大小为5*2
  
B_perp = [eye(4);zeros(1,4)];            %这是那个B垂直

ctrl_N = 12;               %线性插值的个数，这是把这么多飞行状态全模拟了一遍
h_range = linspace(h0-h_lim, h0+h_lim, ctrl_N);      %注意，飞行状态插值应在每个飞行状态附近，即[h0-h_lim,h0+hlim]
v_range = linspace(v0-v_lim, v0+v_lim, ctrl_N);
gama_range = linspace(gama0-gama_lim, gama0+gama_lim, ctrl_N);
alpha_range = linspace(alpha0-alpha_lim, alpha0+alpha_lim, ctrl_N);

%----------------------BEGIN：动力学和微分动力学---------------------------%

%handle function，句柄函数，生成三阶切比雪夫多项式
sin_x_cheby = @(xx) 0.9101*(xx/(pi/3)) - 0.04466*(4*(xx/(pi/3))^3 - 3*(xx/(pi/3)));
cos_x_cheby = @(xx) 0.7441 -0.2499*(2*(xx/(pi/3))^2 -1);

sin_g = @(x) sin_x_cheby(x(4));        
cos_g = @(x) cos_x_cheby(x(4));        %gama
sin_a = @(x) sin_x_cheby(x(5));
cos_a = @(x) cos_x_cheby(x(5));        %alpha

%dynamics f，动力学
S   = @(x) 0.5026;                       %参考面积
hf  = @(x) -x(2)/7110;                   %h fake，中间变量
rou = @(x) 1.225 * (1 + hf + hf^2/2 + hf^3/6 + hf^4/24 + hf^5/120 + hf^6/720 + hf^7/5040 + hf^8/40320 + hf^9/362880 + hf^10/3628800);
q   = @(x) 0.5 * rou * x(3)^2;           %动压
qf  = @(x) 0.5 * rou * x(3);             %q fake，伪动压，已约减速度v
M   = @(x) x(3) / 340;                   %马赫数
m   = @(x) 1600;                         %质量

CL  = @(x) 0.4172 + 19.41*x(5) + 10.17*x(5)^2 - M*(0.1004 + 0.7536*x(5));
L   = @(x) q*CL*S;                       %升力
Lf  = @(x) qf*CL*S;                      %L fake，伪升力，已约减速度v

Cd  = @(x) 0.3042 + 0.02988*CL^2;
D   = @(x) q*Cd*S;                       %阻力

%三阶切比雪夫逼近1/V，区间[1000,2000]
V_division = @(x) 0.0028498843 - 2.984391e-06*x(3) + 1.3615968e-09*x(3)^2 - 2.285664e-13*x(3)^3;

f1 = @(x) -x(3) * cos_a;                 %不是很必要，先留存在这里
f2 = @(x) x(3) * sin_a;
f3 = @(x) -D/m - g*sin_g;
% f4 = L/(m*x(3)) - g*cos_g/x(3);   %(msspoly的除法运算是不被允许的)
% f4 = Lf/m - g*cos_g/x(3);        %Lf已约减速度v，g*cos_g/x(3)还不行
% 如何处理？？先搁置一下，2022/03/27/10:10
f4 = @(x) Lf/m - g*cos_g*V_division;
f5 = @(x) 0;

f_mat = @(x) [f2;f3;f4;f5];               %动力学

%基于微分动力学的df_mat的前四行，控制会影响到后一行，所以拿掉了
df11 = @(x) 0;
df12 = @(x) 0;
df13 = @(x) -cos_a;
df14 = @(x) 0;
df15 = @(x) x(3) * sin_a;
 df1 = @(x) [df11,df12,df13,df14,df15];

df21 = @(x) 0;
df22 = @(x) 0;
df23 = @(x) sin_a;
df24 = @(x) 0;
df25 = @(x) x(3) * cos_a;
 df2 = @(x) [df21,df22,df23,df24,df25];

df31 = @(x) 0;
df32 = @(x) diff(f3,x(2));
df33 = @(x) diff(f3,x(3));
df34 = @(x) diff(f3,x(4));
df35 = @(x) diff(f3,x(5));
 df3 = @(x) [df31,df32,df33,df34,df35];
 
df41 = @(x) 0;
df42 = @(x) diff(f4,x(2));
df43 = @(x) diff(f4,x(3));
df44 = @(x) diff(f4,x(4));
df45 = @(x) diff(f4,x(5));
 df4 = @(x) [df41,df42,df43,df44,df45];
 
df_mat = @(x) [df1;df2;df3;df4;zeros(1,5)];        %微分动力学

%------------------------END：动力学和微分动力学---------------------------%         

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);       %是四维数组，因为有4个状态
eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,2);                     %W只对i,j,k做，因为你的u会影响攻角，所以期望W只是高度、速度、弹道倾角的函数
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);

for i = 1:length(h_range)
    for j = 1:length(v_range)
        for k = 1:length(gama_range)
            for l = 1:length(alpha_range)      %这是小写“L”
                
                x = [randn(1,1);h_range(i);v_range(j);gama_range(k);alpha_range(l)];    %射程随机，飞行状态遍历
                
                W = W_eval(w_poly_fnc(x));
                M = W\eye(n);               %这是在求逆
                Theta = chol(M);            %R = chol(A) 基于矩阵 A 的对角线和上三角形生成上三角矩阵 R，满足方程 R'*R=A。
                Theta_Bw = Theta*Bw(x);
                sigma_ThBw(i,j,k,l) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %为计算偏差上界而准备
                
                L = chol(W);
                
                f = f_mat(x);           %状态微分
                df = df_mat(x);         %状态切向微分
                
                % W是对称的，F也是对称的，本质是G（x）
                F = -W_eval(dw_poly_h_fnc(x)) * f(2) - W_eval(dw_poly_v_fnc(x)) * f(3) - W_eval(dW_gama_fnc(x)) * f(4)...      
                            + df*W + W*df' + 2*lambda*W;  %三个偏微分已经构成了全微分，alpha的全是0，直接不用写
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));     %L是W的分解
                delta_u(i,j,k,l) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %这是控制量
                
                R_CCM = -B_perp'*F*B_perp;  %前后都乘B垂直
                
                eig_CCM(i,j,k,l) = min(eig(R_CCM));
                eig_W(i,j,1) = min(eig(W));
                eig_W(i,j,2) = max(eig(W));
                
            end
        end
    end
end

alpha_w = max(sigma_ThBw(:));       %(:)表示其中所有元素并用列表示
d_bar = alpha_w/lambda;             %这d_bar没有乘w上界，属于是J_CCM也不用再除w上界了，一步到位
disp('d_bar'); 
disp(d_bar);

disp('Control:'); 
disp(max(d_bar*delta_u(:)));        %这是控制上界

disp('W:'); 
disp(min(min(eig_W(:,:,1))));      %小中小
disp(max(max(eig_W(:,:,2))));      %大中大

disp('CCM:'); 
disp(min(eig_CCM(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));
