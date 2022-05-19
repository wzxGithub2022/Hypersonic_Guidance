%% compute_CCM_Hyper线性
% 版本特性与说明：V2 B4
% box_lim应该删不得，这要是删了你的h(x)就没意义了
% 可以考虑停止状态扩展，控制量就用攻角，你的B和B_perp都需要改写
% ddf*x改A*A，以攻角作为控制

% clear all; close all; clc;
% yalmip('clear'); %画图包里的功能，暂且不用
% warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load，导入状态点的信息

load('CCM_Dynamitics_B.mat');
load('CCM_Dynamitics_df.mat');

global state_num state_CCM CCW_upper

% state_num = 1;  %单独跑一个点需要把这一行放出来

%线性化之后信息存在于df_full和B中，state_base/f_ful/ddf_full都没什么用
i = state_num;
B_full = B_mat_value(:,:,i);        % 4*1
df_full = df_mat_value(:,:,i);      % 4*4 double

global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

%% Constants，常数

n = 4;                      %动力学有4个量
ccm_eps = 0.05;             %这个值0.05合适吗？是不是太严苛了

%注意，这是每一点的误差限
h_lim = 100/R0;                  
v_lim = 30/sqrt(R0*g);
gama_lim = pi/180;          %角度的单位是：弧度

%% Uncomment for global optimization of bound，这就是OPT-CCM

% 寻找最小条件数的过程，先1.2倍迭代增大，再二分法搜索
% 优化的终点是什么，在不同的收缩率下，CCM条件数达标时，性能指标是不同的，给一系列lambda，画出对应的性能指标，人工挑选最优的lambda
% euc前期是sqrt(double(w_upper/w_lower))/lambda，等价于（1/lambda^2）*(w_upper/w_lower)，后期是d_bar * diag(W_upper)^1/2
% d_bars前期是sqrt(double(1/w_lower))/lambda，后期是alpha_w/lambda

% lambda_range = linspace(0.1,1,10);              %线性插值，lambda的范围
lambda_range = 0.5;                               %定住了
lambda_range = (1/100)*round(lambda_range*100);   %四舍五入圆整
euc_bounds = NaN(length(lambda_range),1);         %注意，前期优化寻找CCM和后期准确计算时的含义不严格相同
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
        [sos_prob,~,~] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,...
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
        if cond_l > 1000       %为了不让程序进入无限增大条件数的死循环中
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
        condn = (cond_l+cond_u)/2;    %condn只是个中间变量
        fprintf(' cond: %.2f', condn);
        [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,...
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

%     condn_prev = cond_u;            %这行代码是真的疯狂，思路是lambda增大，预设条件数也增大
    cond_bound(ll) = cond_u;        %误差较小，全都相等

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end
% figure(state_num)
% subplot(211),plot(lambda_range,euc_bounds);title('Euc_bound');
% subplot(212),plot(lambda_range,d_bars);title('d_bars');
% pause; 

%% Pick a (lambda,condn)，这是组合拳，挑一套的

lambda = 0.5;          %这是挑出来的
condn = cond_u;        %固定收缩率，条件数自适应
return_metric = 1;

save_name = 'metric_Hyper_vectorized_linear.mat';          %之后的部分只调用了一次CCM_OTP，还用save_name把东西取出来，之前打掉的部分已经完成了优化CCM的功能
%%
[sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,...
                                df_full,condn,lambda,ccm_eps,return_metric,save_name);      %在这一步存好，sos_prob:sum of squares problem

%% Compute aux control bound，辅助控制边界
%这是想看看在算出来的(lambda，condn)情况下，CCM有什么样的表现，需要多大的u才能控得住
%控制和扰动是耦合的，给你一个扰动上界，你还想把它控制在什么样tube内，决定了你的u的情况

load(save_name);
disp('Checking CCM conditions and Computing control bound...');

%----------------------BEGIN：动力学和微分动力学---------------------------%

syms y h v gama           %暂定为归一化的量，写归一化的动力学

%dynamics f，动力学，当你看透了这里的f_mat和后期的f(1...4)的关系，完全不用刻意区分compute和opt的写法区别
f_mat(y, h, v, gama) = df_full * [y; h; v; gama];               %动力学，4*1，含x

%微分动力学的df_mat
df44 = df_full * df_full;
df_mat = df44;            %4*4，是微分系数矩阵，不含delta_x

%------------------------END：动力学和微分动力学---------------------------%

%%
%B只在求控制上界那里用过，可以考虑扔掉alpha_dot，攻角就是控制量，更直接更强有力的控制肯定是更好控的
B = B_full;       % 4*1 的double

%写为0.1倍的气动力
load('CCM_Dynamitics_D.mat')
load('CCM_Dynamitics_Lf.mat')
D_full = D_mat_value(state_num);
Lf_full = Lf_mat_value(state_num);

Bw = [zeros(2,2);
      D_full*0.01, 0;
      0, Lf_full*0.01];
  
B_perp = [eye(2);zeros(2,2)];            %这是那个B垂直
%%
ctrl_N = 6;               %线性插值的个数，这是把这么多飞行状态全模拟了一遍
h_range = linspace(-h_lim, h_lim, ctrl_N);      %注意，飞行状态插值应在每个飞行状态附近，即[h0-h_lim,h0+hlim]
v_range = linspace(-v_lim, v_lim, ctrl_N);
gama_range = linspace(-gama_lim, gama_lim, ctrl_N);

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N);       %是3维数组，因为要遍历3个状态
eig_CCM_min = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_CCM_max = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,2);                     %W只对i做，因为你的攻角会影响v和gama，所以期望W只是高度h的函数
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);
%%
for i = 1:length(h_range)
    fprintf(' i = %d , calculating \n',i);
    for j = 1:length(v_range)
        fprintf(' j = %d , calculating \n',j);
        for k = 1:length(gama_range)
                
                x = [randn(1,1);h_range(i);v_range(j);gama_range(k)];    %射程随机，飞行状态遍历
                
                W = W_eval(w_poly_fnc(x));  %由状态x执行W的生成规则得到W
                M = W\eye(n);               %这是在求逆
                Theta = chol(M);            %R = chol(A) 基于矩阵 A 的对角线和上三角形生成上三角矩阵 R，满足方程 R'*R=A。
                Theta_Bw = Theta*Bw;
                sigma_ThBw(i,j,k) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %为计算偏差上界而准备
                
                L = chol(W);
                
                %和y没关系，写上反正f和df这边乘出来也是0
                f = double( f_mat(0,x(2),x(3),x(4)) );           %状态微分
                df = df_mat;         %状态切向微分
                
                % W是对称的，F也是对称的，本质是G（x），F需要负定，这里的对应关系有过bug，先改为敬
                F = -W_eval(dw_poly_h_fnc(x)) * f(2) + df*W + W*df' + 2*lambda*W;  %B4的全微分就只有一个偏微分了
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));     %L是W的分解
                delta_u(i,j,k) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %这是控制量
                
                R_CCM = -B_perp'*F*B_perp;  %前后都乘B垂直，R_CCM需要正定
                
                eig_CCM_min(i,j,k) = min(eig(R_CCM));
                eig_CCM_max(i,j,k) = max(eig(R_CCM));
                eig_W(i,1) = min(eig(W));
                eig_W(i,2) = max(eig(W));
                
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
disp(min( eig_W(:,1) ));      %小中小
disp(max( eig_W(:,2) ));      %大中大

disp('min eig CCM:'); 
disp(min(eig_CCM_min(:)));

disp('max eig CCM:'); 
disp(max(eig_CCM_max(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));
%% auto test all 需求
state_CCM(1,state_num) = alpha_w/lambda;
state_CCM(2,state_num) = max(d_bar*delta_u(:));
state_CCM(3,state_num) = min(eig_W(:,1));
state_CCM(4,state_num) = max(eig_W(:,2));
state_CCM(5,state_num) = min(eig_CCM_min(:));
state_CCM(6,state_num) = max(eig_CCM_max(:));
state_CCM(7:10,state_num) = d_bar*sqrt(diag(W_upper));
CCW_upper(:,:,state_num) = W_upper;