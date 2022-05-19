%% 简化模型e1
function [solved,w_lower,w_upper] = Opt_Model_e1(n,state_lim,condn,lambda,ccm_eps,return_metric,varargin)

%% Current state load，导入状态点的信息
global R0 g

%% State-space and dynamics，状态空间和动力学

W_scale = diag([0.001;0.001;0.001;0.001;0.001]);
norm_scale = 1e-4;

%states，状态：y,h,v,gama,alpha
x = msspoly('x',5);         %注意：这些x也是归一化后的，算动力学系数时需要反归一化

%pos_def indeterminates，正定的未确定值
dfiv = msspoly('dfiv',5);   %d5，系统共5个状态
dfor = msspoly('dfor',4);   %d4，控制会影响一个状态

%----------------------BEGIN：动力学和微分动力学---------------------------%

%dynamics f，动力学
m  = 1600;        %质量

CL = 0.4 * 180 / pi;
Cd = 0.02988 * CL^2;

%handle function，句柄函数，生成三阶切比雪夫多项式
sin_x = @(x) 0.985400513847416*x - 0.142496853019355*x^3;
cos_x = @(x) 0.999396553656190 - 0.495559134405119*x^2 + 0.036782872656059*x^4;

sin_g = sin_x(x(4));        
cos_g = cos_x(x(4));        %gama

%切比雪夫逼近1/V，对归一化之后的V，区间[3,7]
V_division = 1.101028089323823 - 0.473922941895471*x(3) + 0.099720451270527*x(3)^2 - 0.010265506659430*x(3)^3 + 4.140771839924636e-04*x(3)^4;

f1 = -x(3) * cos_g;
f2 = x(3) * sin_g;
f3 = - Cd/(m*g) * x(5)^2 - sin_g;
f4 = CL/(m*g) * x(5) * V_division - cos_g * V_division;
f5 = 0;

f = [f1;f2;f3;f4;f5];               %动力学，这里的f得和x保持同样的维度

%基于微分动力学的df_mat的前四行，控制会影响到后一行，所以拿掉了
df11 = 0;
df12 = 0;
df13 = -cos_g;
df14 = x(3) * sin_g;
df15 = 0;
 df1 = [df11,df12,df13,df14,df15];

df21 = 0;
df22 = 0;
df23 = sin_g;
df24 = x(3) * cos_g;
df25 = 0;
 df2 = [df21,df22,df23,df24,df25];

df31 = 0;
df32 = 0;
df33 = 0;
df34 = -cos_g;
df35 = - Cd/(m*g) * 2 * x(5);
 df3 = [df31,df32,df33,df34,df35];
 
df41 = 0;
df42 = 0;
df43 = - CL/(m*g) * x(5) * V_division^2 + cos_g * V_division^2;
df44 = sin_g * V_division;
df45 = CL/(m*g) * V_division;
 df4 = [df41,df42,df43,df44,df45];
 
df_perp = [df1;df2;df3;df4];        %微分动力学，这里的df得和B_perp配合

%------------------------END：动力学和微分动力学---------------------------%

B_perp = [eye(4);zeros(1,4)];       %这是那个B垂直

%% Initialize problem，搭建问题框架

prog = spotsosprog;
prog = prog.withIndeterminate(x);           %Indeterminate，待定量
prog = prog.withIndeterminate(dfiv);        %x是状态量
prog = prog.withIndeterminate(dfor);        %dfiv dfor是辅助变量

[prog, w_lower] = prog.newPos(1);           %positive，生成了一个正的值
[prog, w_upper] = prog.newPos(1);

%% Parametrize W，将W参数化

w_states = x(3:4);                      %这是说w作为一个ccm应该是速度v、弹道倾角gama的function，这其实是CCM的条件之一，生成W所需要的状态成为w_states
w_poly = monomials(w_states,0:4);       %三个状态量，应该写出0阶到6阶的单项式，0:6--84，0:4--35

W_list      = cell(length(w_poly),1);
W_perp_list = cell(length(w_poly),1);
W_pc_list   = cell(length(w_poly),1);
W_c_list    = cell(length(w_poly),1);

[prog, W_perp_list{1}] = prog.newSym(4);    %定义了4维的新系统，新的4维对称矩阵symmetric
[prog, W_pc_list{1}] = prog.newFree(4,1);   %newSym是退化的newFree，新的自由变量
[prog, W_c_list{1}] = prog.newSym(1);
W_list{1} = [W_perp_list{1},W_pc_list{1};
             W_pc_list{1}', W_c_list{1}];   %这么大只

W = W_list{1}*w_poly(1);    %做第一个,就是说W和w_poly的每一项都建立关系，而w_poly和w_states的每一项的0到6次单项式建立关系

for i = 2:length(w_poly)
    [prog, W_perp_list{i}] = prog.newSym(4);    %形式上和第一个完全一样
    [prog, W_pc_list{i}] = prog.newFree(4,1);
    [prog, W_c_list{i}] = prog.newSym(1);
    W_list{i} = [W_perp_list{i},W_pc_list{i};
        W_pc_list{i}', W_c_list{i}];
    
    W = W + W_list{i}*w_poly(i);    %并累加
end

W_perp = W(1:4,1:4);                    %去除控制直接影响的
dW_f_perp = diff(W_perp(:),x)*f;        %这个diff是求微分，(:)已经写成一列的形式了，这里的f得和x保持同样的维度
dW_f_perp = reshape(dW_f_perp,4,4);     %吴恩达说过reshape是个好东西

[prog, W_upper] = prog.newSym(n);       %W_upper是5*5msspoly，新的对称矩阵

%% Definiteness conditions，约束条件
%---------------------------------keep box--------------------------------%
v_lim_l = state_lim(1);
v_lim_u = state_lim(2);
v_sum = (v_lim_l + v_lim_u) / 2;
v_minus = (v_lim_u - v_lim_l) / 2;

gama_lim_l = state_lim(3);           %这两个角度的单位都是：弧度
gama_lim_u = state_lim(4);
gama_sum = (gama_lim_l + gama_lim_u) / 2;
gama_minus = (gama_lim_u - gama_lim_l) / 2;

alpha_lim_l = state_lim(5);
alpha_lim_u = state_lim(6);
alpha_sum = (alpha_lim_l + alpha_lim_u) / 2;
alpha_minus = (alpha_lim_u - alpha_lim_l) / 2;
%Lagrange multipliers，拉格朗日乘子
box_lim = [ v_minus^2 - (x(3)-v_sum)^2;
            gama_minus^2 - (x(4)-gama_sum)^2;
            alpha_minus^2 - (x(5)-alpha_sum)^2];       %限定了状态范围
        
%加上跑了也没什么用  
prog = prog.withSOS( box_lim(1) );   
prog = prog.withSOS( box_lim(2) ); 
prog = prog.withSOS( box_lim(3) ); 

%这里是约束W，而W只和w_states有关
l_order = 4;    %与控制无关的状态2阶，dfiv需要2阶
[pos_def_mon, pos_def_mat] = monomials([w_states;dfiv],0:l_order);      %monomials单项式，mat记录了生成规则
pos_def_keep = find( sum( pos_def_mat(:, length(w_states)+1:end) ,2 ) == 2 );   %保留dfiv的2阶
pos_def_mon = pos_def_mon(pos_def_keep);            
[prog, Ll] = prog.newSOSPoly(pos_def_mon,2);
[prog, Lu] = prog.newSOSPoly(pos_def_mon,2);        %L1 Lu都是2*1的msspoly，小L是低值，u是高值，生成了新的SOS的poly

%这里是约束G(x)/R_CCM，与影响系统的状态量都有关
lc_order = 5;   %包含控制相关的状态3阶，dfor需要2阶
l_ccm_states = [x(3);x(4);x(5)];
[ccm_def_mon, ccm_def_mat] = monomials([l_ccm_states;dfor],0:lc_order);          
ccm_def_keep = find(sum(ccm_def_mat(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dfor
ccm_def_mon = ccm_def_mon(ccm_def_keep);
[prog, Lc] = prog.newSOSPoly(ccm_def_mon,3);

%W uniform bounds，W的范围约束
prog = prog.withPos(w_lower-1);                     %w_lower > 1，正条件约束
prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper，正定条件约束

%Condition bound，条件数约束
prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn，正条件约束

%W pos def，positive definite，正定约束，用SOS的方式来约束，在box范围内W都满足条件
prog = prog.withSOS( (dfiv'*W*dfiv - w_lower*(dfiv'*dfiv)) - (Ll'*box_lim(1:2)) );      %约束 W > w_lower*eye5
prog = prog.withSOS( dfiv'*(W_upper - W)*dfiv - (Lu'*box_lim(1:2)));                    %约束 W_upper > W

%CCM condition，CCM条件，保证收缩
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %本质是G（x）那个条件，前后都乘个B垂直，W就变成了W_perp
prog = prog.withSOS( (dfor'*R_CCM*dfor - ccm_eps*(dfor'*dfor)) - (Lc'*box_lim) );   %约束 R_CCM > epsilion*eye4
%---------------------------------keep box--------------------------------%

%--------------------------------throw box--------------------------------%
% %W uniform bounds，W的范围约束
% prog = prog.withPos(w_lower-1);                     %w_lower > 1，正条件约束
% prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper，正定条件约束
% %Condition bound，条件数约束
% prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn，正条件约束
% %W pos def，positive definite，正定约束，用SOS的方式来约束，在box范围内W都满足条件
% prog = prog.withSOS( (dfiv'*W*dfiv - w_lower*(dfiv'*dfiv)) );      %约束 W > w_lower*eye5
% prog = prog.withSOS( dfiv'*(W_upper - W)*dfiv );                    %约束 W_upper > W
% %CCM condition，CCM条件，保证收缩
% R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %本质是G（x）那个条件，前后都乘个B垂直，W就变成了W_perp
% prog = prog.withSOS( (dfor'*R_CCM*dfor - ccm_eps*(dfor'*dfor)) );   %约束 R_CCM > epsilion*eye4
%--------------------------------throw box--------------------------------%

%约束完成，构建SDP半定规划问题
options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint，常量约束，压缩解的范围，使问题更易求解
free_vars = [prog.coneVar(2:end); prog.freeVar];   
len = length(free_vars);
[prog, a] = prog.newPos(len);       %a和free_vars具有相同的维度，正条件约束
prog = prog.withPos(-free_vars + a);    
prog = prog.withPos(free_vars + a);

%try-catch-end，执行语句并捕获产生的错误
    SOS_soln = prog.minimize(trace(W_scale*W_upper) + norm_scale*sum(a), @spot_mosek, options);     %调用mosek，直接一个minimize优化结束了
    %不归一化的话，这里的优化应该跑不出结果，算力障碍本质体现为问题设计不妥
try
    %SOS_soln就是prog
    solved = ~SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');    %比较字符串，当比较成功时，solved会变成0
    w_lower = double(SOS_soln.eval(w_lower));
    w_upper = double(SOS_soln.eval(w_upper));
catch
    %failed，这意味着CCM失败，不能解决问题
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    return;
end

%% Parse，解析，把结果提取出来

W_upper_mat = zeros(n);     %属实没有用得上，在3D-Quadrotor用得到

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');

        W_sol = zeros(n,n,length(w_poly));  %从SOS_soln.eval(W_list{i})里取值并筛掉小量
        NNZ_list = zeros(length(w_poly),1); 
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-9); %clean Remove terms with small coefficients
            if sum(sum(abs(W_sol(:,:,i)))) > 0      %绝对值大于零，证明该W_sol(:,:,i)算出东西了，保留
                NNZ_list(i) = 1;                    %non_negative，保留的一个标记
            end
        end
        w_poly = w_poly(find(NNZ_list));        %find返回的是索引，w_poly重获新生
        W_sol = W_sol(:,:,find(NNZ_list));      %W_sol也重获新生，好像也没排除什么
        
        fprintf('%d non-zero monomials\n',length(w_poly)); %统计非零（极小值tolerance）单项
        
        dw_poly_v = diff(w_poly,x(3));
        dw_poly_gama = diff(w_poly,x(4));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-9);

        %% Create monomial functions，创造句柄函数
        w_poly_fnc       = mss2fnc(w_poly,    x,randn(length(x),2));  %computes all monomials
        dw_poly_v_fnc    = mss2fnc(dw_poly_v, x,randn(length(x),2));  %理论上dw_h + dw_v + dw_gama = w
        dw_poly_gama_fnc = mss2fnc(dw_poly_gama,x,randn(length(x),2));  %就是说三个偏导加起来就是全微分

        %% Put together，累加句柄
        W_exec = 'W_eval = @(ml)';  %W_eval的执行对象(m1)是以上三个求偏导
        
        for i = 1:length(w_poly)
            if i<length(w_poly)
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d) +',i,i));      %水平串联字符串
            else
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d);',i,i));
            end
        end

        %% Execute，算好了保存输出
        eval(W_exec);   %computes the sum above，执行文本中的 MATLAB 表达式，文本转功能
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_v_fnc','dw_poly_gama_fnc','W_upper'); %都存进去，以供调用

    end
end
end
