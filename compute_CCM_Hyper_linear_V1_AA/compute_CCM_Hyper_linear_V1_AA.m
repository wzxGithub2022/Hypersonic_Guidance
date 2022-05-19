%% compute_CCM_Hyper线性
% 版本特性与说明：V1 AA
% 主程序里完全没有用到spotless的东西，他们只是优化用的工具
% 问题还出在opt里
% box_lim应该删不得，这要是删了你的h(x)就没意义了
% 可以考虑停止状态扩展，控制量就用攻角，你的B和B_perp都需要改写
% ddf*x改A*A，仍以攻角变化率作为控制
clear all; close all; clc;
% yalmip('clear'); %画图包里的功能，暂且不用
warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load，导入状态点的信息

load('Trajectory_normalization.mat');
load('CCM_Dynamitics_df.mat');

%线性化之后信息存在于df_full和ddf_full中，state_base和f_full都没什么用
i = 1;
df_full = df_mat_value(:,:,i);      % 5*5 double

global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

%% Constants，常数

n = 5;                      %动力学有5个量
ccm_eps = 0.001;             %这个值0.05合适吗？是不是太严苛了

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

lambda_range = linspace(0.1,2,10);              %线性插值，lambda的范围
lambda_range = (1/100)*round(lambda_range*100);   %四舍五入为整数
euc_bounds = NaN(length(lambda_range),1);         %注意，前期优化寻找CCM和后期准确计算时的含义不严格相同，
d_bars = NaN(length(lambda_range),1);             %bar是上面一横
cond_bound = NaN(length(lambda_range),1);         %条件数

eps = 1;
condn_prev = 5;                                  %预设条件数
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
        [sos_prob,~,~] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,alpha_lim,...
                                df_full,cond_u,lambda,ccm_eps,return_metric);
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
        if cond_l > 1000
            break;
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
        [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,alpha_lim,...
                                df_full,condn,lambda,ccm_eps,return_metric);
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

%     condn_prev = cond_u;      %这行的逻辑先停掉         
    cond_bound(ll) = cond_u;    %误差较小，全都相等

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end

pause; 

% %% Pick a (lambda,condn)，这是组合拳，挑一套的
% 
% lambda = 0.7;          %这是挑出来的
% condn = 10;
% return_metric = 1;
% 
% save_name = 'metric_Hyper_vectorized_linear.mat';          %之后的部分只调用了一次CCM_OTP，还用save_name把东西取出来，之前打掉的部分已经完成了优化CCM的功能
% %%
% [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,alpha_lim,...
%                                 df_full,condn,lambda,ccm_eps,return_metric,save_name);      %在这一步存好，sos_prob:sum of squares problem
% 
% %% Compute aux control bound，辅助控制边界
% %这是想看看在算出来的(lambda，condn)情况下，CCM有什么样的表现，需要多大的u才能控得住
% %控制和扰动是耦合的，给你一个扰动上界，你还想把它控制在什么样tube内，决定了你的u的情况
% 
% load(save_name);
% disp('Checking CCM conditions and Computing control bound...');
% 
% %B只在求控制上界那里用过，可以考虑扔掉alpha_dot，攻角就是控制量，更直接更强有力的控制肯定是更好控的
% B = [zeros(4,1);1];       % 5*1 的double
%  
% Bw = @(x)[zeros(2,2);     % 这两个扰动的意思是：升力和阻力误差为0.01倍的重力
%            0.01/1600,0;
%      0,0.01/(1600*2000);
%                 0,0];      %这样的写法直接就是小function啊，气动力扰动对系统的输入，大小为5*2
%   
% B_perp = [eye(4);zeros(1,4)];            %这是那个B垂直
% 
% ctrl_N = 6;               %线性插值的个数，这是把这么多飞行状态全模拟了一遍
% h_range = linspace(-h_lim, h_lim, ctrl_N);      %注意，飞行状态插值应在每个飞行状态附近，即[h0-h_lim,h0+hlim]
% v_range = linspace(-v_lim, v_lim, ctrl_N);
% gama_range = linspace(-gama_lim, gama_lim, ctrl_N);
% alpha_range = linspace(-alpha_lim, alpha_lim, ctrl_N);
% 
% %----------------------BEGIN：动力学和微分动力学---------------------------%
% 
% syms h v gama alpha           %暂定为归一化的量，写归一化的动力学
% 
% % %dynamics f，动力学
% f_mat(h, v, gama, alpha) = df_full(2:5,2:5) * [h; v; gama; alpha];               %动力学，4*1，含x
% 
% %基于微分动力学的df_mat的前四行，控制会影响到后一行，所以拿掉了
% df55 = df_full(:,:) * df_full(:,:) * [0; h; v; gama; alpha];
% df_mat(h, v, gama, alpha) = df55;            %5*5，含x
% 
% %------------------------END：动力学和微分动力学---------------------------%
% 
% delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);       %是四维数组，因为有4个状态
% eig_CCM_min = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
% eig_CCM_max = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
% eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,2);                     %W只对i,j,k做，因为你的u会影响攻角，所以期望W只是高度、速度、弹道倾角的函数
% sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);
% %%
% for i = 1:length(h_range)
%     fprintf(' i = %d , calculating \n',i);
%     for j = 1:length(v_range)
%         fprintf(' j = %d , calculating \n',j);
%         for k = 1:length(gama_range)
%             for l = 1:length(alpha_range)      %这是小写“L”
%                 
%                 x = [randn(1,1);h_range(i);v_range(j);gama_range(k);alpha_range(l)];    %射程随机，飞行状态遍历
%                 
%                 W = W_eval(w_poly_fnc(x));
%                 M = W\eye(n);               %这是在求逆
%                 Theta = chol(M);            %R = chol(A) 基于矩阵 A 的对角线和上三角形生成上三角矩阵 R，满足方程 R'*R=A。
%                 Theta_Bw = Theta*Bw(x);
%                 sigma_ThBw(i,j,k,l) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %为计算偏差上界而准备
%                 
%                 L_nor = chol(W);
%                 
%                 f = double( f_mat(x(2),x(3),x(4),x(5)) );           %状态微分
%                 df = double( df_mat(x(2),x(3),x(4),x(5)) );         %状态切向微分
%                 
%                 % W是对称的，F也是对称的，本质是G（x），F需要负定，先改为敬
%                 F = -W_eval(dw_poly_h_fnc(x)) * f(1) - W_eval(dw_poly_v_fnc(x)) * f(2) - W_eval(dw_poly_gama_fnc(x)) * f(3)...      
%                             + df*W + W*df' + 2*lambda*W;  %三个偏微分已经构成了全微分，alpha的全是0，直接不用写
%                 
%                 delta_u_den = eig((inv(L_nor))'*(B*B')*inv(L_nor));     %L是W的分解
%                 delta_u(i,j,k,l) = 0.5*max(eig((inv(L_nor))'*F*inv(L_nor)))/...
%                     sqrt(min( delta_u_den(delta_u_den>0) ));    %这是控制量
%                 
%                 R_CCM = -B_perp'*F*B_perp;  %前后都乘B垂直，R_CCM需要正定
%                 
%                 eig_CCM_min(i,j,k,l) = min(eig(R_CCM));
%                 eig_CCM_max(i,j,k,l) = max(eig(R_CCM));
%                 eig_W(i,j,k,1) = min(eig(W));
%                 eig_W(i,j,k,2) = max(eig(W));
%                 
%             end
%         end
%     end
% end
% %%
% alpha_w = max(sigma_ThBw(:));       %(:)表示其中所有元素并用列表示
% d_bar = alpha_w/lambda;             %这d_bar没有乘w上界，属于是J_CCM也不用再除w上界了，一步到位
% disp('d_bar'); 
% disp(d_bar);
% 
% disp('Control:'); 
% disp(max(d_bar*delta_u(:)));        %这是控制上界
% 
% disp('W:'); 
% disp(min(min(min(eig_W(:,:,:,1)))));      %小中小
% disp(min(max(max(eig_W(:,:,:,2)))));      %大中大
% 
% disp('min eig CCM:'); 
% disp(min(eig_CCM_min(:)));
% 
% disp('max eig CCM:'); 
% disp(max(eig_CCM_max(:)));
% 
% disp('euc_bounds');
% disp(d_bar*sqrt(diag(W_upper)));
