%% compute_CCM_Hyper
% 版本特性与说明：V1 CAD control_alpha_dot 
% 离谱的alpha bug已修改
% 当前以攻角变化率作为控制量，思路是科学的，但算力障碍且给(0.7,10)时G(x)不满足收敛条件，理论上可以抢救
% 考虑停止状态扩展，控制量用攻角，B和B_perp都需要改写，生成V2 CA
% 可以调整的参数：ccm_eps, normscale

clear all; close all; clc;
% yalmip('clear'); %画图包里的功能，暂且不用
warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load，导入状态点的信息

load('Trajectory_normalization.mat');

% global state_num state_CCM

i_state = 1;
h1 = Trajectory_normalization(i_state,3);
v1 = Trajectory_normalization(i_state,4);
gama1 = Trajectory_normalization(i_state,5);
alpha1 = Trajectory_normalization(i_state,6);
state_base = [h1,v1,gama1,alpha1];

global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

%% Constants，常数

n = 5;                      %动力学有5个量
ccm_eps = 0.01;

%注意，这是每一点的误差限
h_lim = 100/R0;                  
v_lim = 30/sqrt(R0*g);
gama_lim = pi/180;           %这两个角度的单位都是：弧度
alpha_lim = pi/360;

%% Uncomment for global optimization of bound，这就是OPT-CCM

%寻找最小条件数的过程，先1.2倍迭代增大，再二分法搜索
%优化的终点是什么，在不同的收缩率下，CCM条件数达标时，性能指标是不同的，给一系列lambda，画出对应的性能指标，人工挑选最优的lambda
%euc前期是sqrt(double(w_upper/w_lower))/lambda，等价于（1/lambda^2）*(w_upper/w_lower)，后期是d_bar * diag(W_upper)^1/2
%d_bars前期是sqrt(double(1/w_lower))/lambda，后期是alpha_w/lambda

% lambda_range = linspace(0.7,0.95,5);              %线性插值，lambda的范围
% lambda_range = (1/100)*round(lambda_range*100);   %四舍五入为整数
% euc_bounds = NaN(length(lambda_range),1);         %注意，前期优化寻找CCM和后期准确计算时的含义不严格相同，
% d_bars = NaN(length(lambda_range),1);             %bar是上面一横
% cond_bound = NaN(length(lambda_range),1);         %条件数
% 
% eps = 1;
% condn_prev = 10;                                  %预设条件数
% return_metric = 0;
% 
% for ll = 1:length(lambda_range)       %对每一个lambda操作，是小写的“LL”,美观，接受了
%     
%     lambda = lambda_range(ll);
%     fprintf('**********\n');
%     fprintf('lambda: %f\n', lambda);
%     solved = 0;
%     
%     %Determine upper bound，决定上界
%     cond_l = condn_prev;              
%     cond_u = 1.2*condn_prev;
% %------------------------若代码一直在循环就是卡在这个区间，sos问题一直得不到解决-----------------------%
%     while (~solved) 
%         fprintf(' cond_u: %.2f: ', cond_u);
%         [sos_prob,~,~] = CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
%                                 cond_u,lambda,ccm_eps,return_metric);
%         if (sos_prob == 0)
%             solved = 1;
%             fprintf('feasible \n');
%         else
%             %shift up condition number range
%             %解决不了的话就一点点调大condition number，只调大不一定解决问题
%             fprintf('你看这解不出来了吧\n');
%             cond_l = cond_u;
%             cond_u = 1.2*cond_u;
%         end
%     end
% %---------------------------------------------------------------------------------------------------%
%     if (solved)           %直到解决了，把euc_bounds(ll)定下来
%         euc_bounds(ll) = sqrt(cond_u)/lambda;       %初步算出来一个值，此时还无法计算d_bars
%         fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
%     else
%         continue;
%     end
%     
%     %Now do bisection search，做二分搜索
%     while(cond_u - cond_l >= eps)     %误差较大，细分优化
%         condn = (cond_l+cond_u)/2;
%         fprintf(' cond: %.2f', condn);
%         [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
%                                 condn,lambda,ccm_eps,return_metric);
%         if (sos_prob == 0)
%             fprintf(' feasible\n');       %这就可行了，向左再取小一点，向左动的话，euc和d就还能再优化一些      
%             euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;   %euc等价于（1/lambda^2）*(w_upper/w_lower)
%             d_bars(ll) = sqrt(double(1/w_lower))/lambda;             %这里才首次计算d_bars
%             cond_u = condn;
%         else
%             fprintf(' infeasible\n');     %不可行的话，向右放大一点
%             cond_l = condn;
%         end
%     end
% 
%     condn_prev = cond_u;              %误差较小，全都相等
%     cond_bound(ll) = cond_u;
% 
%     disp('Euc_bound:'); disp(euc_bounds(ll));
%     disp('d_bar:'); disp(d_bars(ll));
%     fprintf('**********\n');
%     
% end
% 
% pause; 

%% Pick a (lambda,condn)，这是组合拳，挑一套的

lambda = 0.7;          %这是挑出来的
condn = 10;
return_metric = 1;

save_name = 'metric_Hyper_vectorized.mat';          %之后的部分只调用了一次CCM_OTP，还用save_name把东西取出来，之前打掉的部分已经完成了优化CCM的功能
%%
[sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                condn,lambda,ccm_eps,return_metric,save_name);      %在这一步存好，sos_prob:sum of squares problem

%% Compute aux control bound，辅助控制边界
%这是想看看在算出来的(lambda，condn)情况下，CCM有什么样的表现，需要多大的u才能控得住
%控制和扰动是耦合的，给你一个扰动上界，你还想把它控制在什么样tube内，决定了你的u的情况
load(save_name);
disp('Checking CCM conditions and Computing control bound...');
%%
%----------------------BEGIN：动力学和微分动力学---------------------------%
% 按照PVTOL，这边是可以用精确的，岔开搭配，msspoly是被迫多项式的
syms h v gama alpha           %暂定为归一化的量，写归一化的动力学

%dynamics f，动力学
S   = 0.5026;                       %参考面积
% rou = 1.225 * exp(-x(2)/7110);    %密度rou，需要多项式逼近
hf  = -h*R0/7110;                    %h fake，中间变量，exp_cheby_bug已改标注
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
f5 = 0;

f_mat(h, v, gama, alpha) = [f2;f3;f4;f5];               %动力学，f存在的意义是做 偏W/偏f 那一项

%基于微分动力学的df_mat的前四行，控制会影响到后一行，所以拿掉了
df1_y = 0;
df1_h = 0;
df1_v = -cos_g;
df1_gama = v * sin_g;
df1_alpha = 0;

df1 = [df1_y, df1_h, df1_v, df1_gama, df1_alpha];

df2_y = 0;
df2_h = 0;
df2_v = sin_g;
df2_gama = v * cos_g;
df2_alpha = 0;

df2 = [df2_y, df2_h, df2_v, df2_gama, df2_alpha];

df3_y = 0;
df3_h = diff(f3,h);
df3_v = diff(f3,v);
df3_gama = -cos_g;
df3_alpha = diff(f3,alpha);

df3 = [df3_y, df3_h, df3_v, df3_gama, df3_alpha];

df4_y = 0;
df4_h = diff(f4,h);
df4_v = diff(f4,v);
df4_gama = sin_g*V_division;
df4_alpha = diff(f4,alpha);

df4 = [df4_y, df4_h, df4_v, df4_gama, df4_alpha];

df5 = [0,0,0,0,0];

df_mat(h, v, gama, alpha) = [df1;df2;df3;df4;df5];      %df都用得上，5*5没毛病

%------------------------END：动力学和微分动力学---------------------------%
%%
B = [zeros(4,1);1];       % 5*1 的double
 
% Bw = @(x)[zeros(2,2);     % 这两个扰动的意思是：升力和阻力误差为0.01倍的重力
%            0.01/1600,0;
%      0,0.01/(1600*2000);
%                 0,0];      %这样的写法直接就是小function啊，气动力扰动对系统的输入，大小为5*2
            
Bw_symfnc(h, v, gama, alpha) = [zeros(2,2);         %气动系数上的扰动引起的Tube就大得离谱
                                0.01*D_nor, 0;
                                0, 0.01*Lf_nor;
                                0,0];     %扰动是0.1倍的气动系数误差
% Bw = double( Bw_symfnc(h1, v1, gama1, alpha1) );
  
B_perp = [eye(4);zeros(1,4)];            %这是那个B垂直
%%
ctrl_N = 6;               %线性插值的个数，这是把这么多飞行状态全模拟了一遍
h_range = linspace(h1-h_lim, h1+h_lim, ctrl_N);      %注意，飞行状态插值应在每个飞行状态附近，即[h0-h_lim,h0+hlim]
v_range = linspace(v1-v_lim, v1+v_lim, ctrl_N);
gama_range = linspace(gama1-gama_lim, gama1+gama_lim, ctrl_N);
alpha_range = linspace(alpha1-alpha_lim, alpha1+alpha_lim, ctrl_N);

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);       %是四维数组，因为有4个状态
eig_CCM_min = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_CCM_max = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,2);                     %W只对i,j,k做，因为你的u会影响攻角，所以期望W只是高度、速度、弹道倾角的函数
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);
%%
for i = 1:length(h_range)
    fprintf(' i = %d , calculating \n',i);
    for j = 1:length(v_range)
        fprintf(' j = %d , calculating \n',j);
        for k = 1:length(gama_range)
            for l = 1:length(alpha_range)      %这是小写“L”
                
                x = [randn(1,1);h_range(i);v_range(j);gama_range(k);alpha_range(l)];    %射程随机，飞行状态遍历
                
                Bw = double( Bw_symfnc(x(2),x(3),x(4),x(5)) );
                
                W = W_eval(w_poly_fnc(x));
                M = W\eye(n);               %这是在求逆
                Theta = chol(M);            %R = chol(A) 基于矩阵 A 的对角线和上三角形生成上三角矩阵 R，满足方程 R'*R=A。
                Theta_Bw = Theta*Bw;
                sigma_ThBw(i,j,k,l) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %为计算偏差上界而准备
                
                L = chol(W);        %注意，L_nor和归一化升力命名冲突，这是你shift+enter的结果，原本就是L
                
                f = double( f_mat(x(2),x(3),x(4),x(5)) );           %状态微分
                df = double( df_mat(x(2),x(3),x(4),x(5)) );         %状态切向微分
                
                % W是对称的，F也是对称的，本质是G（x），F需要负定，这里的对应关系好像有问题？
                F = -W_eval(dw_poly_h_fnc(x)) * f(1) - W_eval(dw_poly_v_fnc(x)) * f(2) - W_eval(dw_poly_gama_fnc(x)) * f(3)...      
                            + df*W + W*df' + 2*lambda*W;  %三个偏微分已经构成了全微分，alpha的全是0，直接不用写
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));     %L是W的分解
                delta_u(i,j,k,l) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %这是控制量
                
                R_CCM = -B_perp'*F*B_perp;  %前后都乘B垂直，R_CCM需要正定
                
                eig_CCM_min(i,j,k,l) = min(eig(R_CCM));     %若他为负值，则不满足收缩条件
                eig_CCM_max(i,j,k,l) = max(eig(R_CCM));
                eig_W(i,j,k,1) = min(eig(W));
                eig_W(i,j,k,2) = max(eig(W));
                
            end
        end
    end
end
%%
alpha_w = max(sigma_ThBw(:));       %(:)表示其中所有元素并用列表示
d_bar = alpha_w/lambda;             %这d_bar没有乘w上界，属于是J_CCM也不用再除w上界了，一步到位
disp('d_bar'); 
disp(d_bar);

disp('Control:'); 
disp(max(d_bar*delta_u(:)));        %这是控制上界

disp('W:'); 
disp(min(min(min(eig_W(:,:,:,1)))));      %小中小
disp(max(max(max(eig_W(:,:,:,2)))));      %大中大

disp('min eig CCM:'); 
disp(min(eig_CCM_min(:)));

disp('max eig CCM:'); 
disp(max(eig_CCM_max(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));

%%
% state_CCM(1,state_num) = alpha_w/lambda;
% state_CCM(2,state_num) = max(d_bar*delta_u(:));
% state_CCM(3,state_num) = min(min(min(eig_W(:,:,:,1))));
% state_CCM(4,state_num) = max(max(max(eig_W(:,:,:,2))));
% state_CCM(5,state_num) = min(eig_CCM_min(:));
% state_CCM(6,state_num) = max(eig_CCM_max(:));
% state_CCM(7:11,state_num) = d_bar*sqrt(diag(W_upper));
%% CAD版本的问题就是CCM收敛条件不行，下面是0.01气动跑出来的结果，2022/4/20
% d_bar
%     0.0191
% 
% Control:
%     0.1288
% 
% W:
%     0.9913
% 
%     1.0474
% 
% min eig CCM:
%    -8.4467
% 
% max eig CCM:
%     5.9149
% 
% euc_bounds
%     0.0190
%     0.0195
%     0.0191
%     0.0190
%     0.0190