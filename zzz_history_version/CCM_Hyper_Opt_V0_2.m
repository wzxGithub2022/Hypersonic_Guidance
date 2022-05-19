%% 优化得到CCM的function
% solved是一个表示SOS问题是否解决的flag
% 用的框架全是spotless，只是调用了Mosek求解器而已
function [solved,w_lower,w_upper] = ...
    CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
    condn,lambda,ccm_eps,return_metric,varargin)
%% Current state load，导入状态点的信息
h1 = state_base(1);
v1 = state_base(2);
gama1 = state_base(3);
alpha1 = state_base(4);

global R0 g
R0 = 10*10^3;                    %R0单位：m
g = 9.81;

h0 = h1*R0;                      %反归一化示例
v0 = v1*sqrt(R0*g);

%% State-space and dynamics，状态空间和动力学
W_scale = diag([0.001;0.001;0.001;0.001;0.001]);
norm_scale = 1e-4;

%handle function，句柄函数，生成三阶切比雪夫多项式
sin_x = @(x) 0.985400513847416*x - 0.142496853019355*x^3;
cos_x = @(x) 0.999396553656190 - 0.495559134405119*x^2 + 0.036782872656059*x^4;

%states，状态：y,h,v,gama,alpha
x = msspoly('x',5);         %注意：这些x也是归一化后的，算动力学系数时需要反归一化

%pos_def indeterminates，正定的未确定值
dfiv = msspoly('dfiv',5);   %d5，系统共5个状态
dfor = msspoly('dfor',4);   %d4，控制会影响一个状态

sin_g = sin_x(x(4));        
cos_g = cos_x(x(4));        %gama
sin_a = sin_x(x(5));
cos_a = cos_x(x(5));        %alpha

%dynamics f，动力学
S   = 0.5026;                       %参考面积
% rou = 1.225 * exp(-x(2)/7110);    %密度rou，(msspoly不能计算指数，泰勒展开逼近)
hf  = x(2)*R0/7110;                 %h fake，中间变量
rou = 1.225 * (1.124508748184077 + 0.459262515874491*hf + 3.155718757366057*hf^2 + 1.426742501850045*hf^3 + 0.374836249394692*hf^4);
q   = 0.5 * rou * (x(3)*sqrt(R0*g))^2;           %动压
qf  = 0.5 * rou * x(3)*sqrt(R0*g);               %q fake，伪动压，已约减速度v
M   = x(3)*sqrt(R0*g) / 340;                     %马赫数
m   = 1600;                                      %质量

CL  = 0.4172 + 19.41*x(5) + 10.17*x(5)^2 - M*(0.1004 + 0.7536*x(5));
L_nor = q*CL*S / (m*g);                        %升力
Lf_nor = qf*CL*S / (m*g);                      %L fake，伪升力，已约减速度v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %阻力

%切比雪夫逼近1/V，对归一化之后的V，区间[3,7]
V_division = 1.101028089323823 - 0.473922941895471*x(3) + 0.099720451270527*x(3)^2 - 0.010265506659430*x(3)^3 + 4.140771839924636e-04*x(3)^4;

f1 = -x(3) * cos_a;
f2 = x(3) * sin_a;
f3 = -D_nor - sin_g;
f4 = Lf_nor - cos_g*V_division;
f5 = 0;

f = [f1;f2;f3;f4;f5];               %动力学

%基于微分动力学的df_mat的前四行，控制会影响到后一行，所以拿掉了
df11 = 0;
df12 = 0;
df13 = -cos_a;
df14 = 0;
df15 = x(3) * sin_a;
 df1 = [df11,df12,df13,df14,df15];

df21 = 0;
df22 = 0;
df23 = sin_a;
df24 = 0;
df25 = x(3) * cos_a;
 df2 = [df21,df22,df23,df24,df25];

df31 = 0;
df32 = diff(f3,x(2));
df33 = diff(f3,x(3));
df34 = diff(f3,x(4));
df35 = diff(f3,x(5));
 df3 = [df31,df32,df33,df34,df35];
 
df41 = 0;
df42 = diff(f4,x(2));
df43 = diff(f4,x(3));
df44 = diff(f4,x(4));
df45 = diff(f4,x(5));
 df4 = [df41,df42,df43,df44,df45];
 
df_perp = [df1;df2;df3;df4];        %微分动力学

B_perp = [eye(4);zeros(1,4)];       %这是那个B垂直

%% Initialize problem，搭建问题框架

prog = spotsosprog;
prog = prog.withIndeterminate(x);           %Indeterminate，待定量
prog = prog.withIndeterminate(dfiv);
prog = prog.withIndeterminate(dfor);

[prog, w_lower] = prog.newPos(1);           %positive，生成了一个正的值
[prog, w_upper] = prog.newPos(1);

%% Parametrize W，将W参数化

w_states = x(2:4);                      %这是说w作为一个ccm应该是高度h、速度v、弹道倾角gama的function，这其实是CCM的条件之一，生成W所需要的状态成为w_states
w_poly = monomials(w_states,0:6);       %三个状态量，应该写出0阶到6阶的单项式

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
dW_f_perp = diff(W_perp(:),x)*f;        %这个diff是求微分，(:)已经写成一列的形式了
dW_f_perp = reshape(dW_f_perp,4,4);     %吴恩达说过reshape是个好东西

[prog, W_upper] = prog.newSym(n);       %W_upper是5*5msspoly，新的对称矩阵

%% Definiteness conditions，约束条件
% 这部分的语言规范应该很大程度上参考mosek

%注意：如果是逐点计算CCM，是否应该在这里的box引入每个状态平衡点，即h_lim^2 - (x(2)-x_standard)^2
%Lagrange multipliers，拉格朗日乘子
box_lim = [ h_lim^2 - (x(2)-h1)^2;
            v_lim^2 - (x(3)-v1)^2;
            gama_lim^2 - (x(4)-gama1)^2;
            alpha_lim^2 - (x(5)-alpha1)^2];       %限定了状态范围

l_order = 5;    %与控制无关的状态三阶，dfiv需要2阶
[pos_def_mon, pos_def_mat] = monomials([w_states;dfiv],0:l_order);      %monomials单项式，mat记录了生成规则
pos_def_keep = find( sum( pos_def_mat(:, length(w_states)+1:end) ,2 ) == 2 );   %保留dfiv的2阶
pos_def_mon = pos_def_mon(pos_def_keep);            
[prog, Ll] = prog.newSOSPoly(pos_def_mon,3);
[prog, Lu] = prog.newSOSPoly(pos_def_mon,3);        %L1 Lu都是2*1的msspoly，小L是低值，u是高值，生成了新的SOS的poly

lc_order = 6;   %包含控制相关的状态三阶，dfiv需要2阶
l_ccm_states = [x(2);x(3);x(4);x(5)];
[ccm_def_mon, ccm_def_mat] = monomials([l_ccm_states;dfor],0:lc_order);          
ccm_def_keep = find(sum(ccm_def_mat(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dfor
ccm_def_mon = ccm_def_mon(ccm_def_keep);
[prog, Lc] = prog.newSOSPoly(ccm_def_mon,4);

%W uniform bounds，W的范围约束
prog = prog.withPos(w_lower-1);                     %w_lower > 1，正条件约束
prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper，正定条件约束

%Condition bound，条件数约束
prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn，正条件约束

%W pos def，positive definite，正定约束
prog = prog.withSOS( (dfiv'*W*dfiv - w_lower*(dfiv'*dfiv)) - (Ll'*box_lim(1:3)) );      %约束为SOS
prog = prog.withSOS( dfiv'*(W_upper - W)*dfiv - (Lu'*box_lim(1:3)));                    %约束为SOS

%CCM condition，CCM条件，保证收缩
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %本质是G（x）那个条件，前后都乘个B垂直，W就变成了W_perp
prog = prog.withSOS( (dfor'*R_CCM*dfor - ccm_eps*(dfor'*dfor)) - (Lc'*box_lim) );   %约束为SOS

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

        W_sol = zeros(n,n,length(w_poly));  
        NNZ_list = zeros(length(w_poly),1); 
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-3); %clean Remove terms with small coefficients
            if sum(sum(abs(W_sol(:,:,i)))) > 0      %绝对值大于零，证明该W_sol(:,:,i)算出东西了，保留
                NNZ_list(i) = 1;                    %non_negative，保留的一个标记
            end
        end
        w_poly = w_poly(find(NNZ_list));        %find返回的是索引，w_poly重获新生
        W_sol = W_sol(:,:,find(NNZ_list));      %W_sol也重获新生，好像也没排除什么
        
        fprintf('%d non-zero monomials\n',length(w_poly)); %统计非零（极小值tolerance）单项
        
        dw_poly_h = diff(w_poly,x(2));   %看这样子就是求偏导
        dw_poly_v = diff(w_poly,x(3));
        dw_poly_gama = diff(w_poly,x(4));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-3);

        %% Create monomial functions，创造句柄函数
        w_poly_fnc       = mss2fnc(w_poly,    x,randn(length(x),2));  %computes all monomials
        dw_poly_h_fnc    = mss2fnc(dw_poly_h, x,randn(length(x),2));
        dw_poly_v_fnc    = mss2fnc(dw_poly_v, x,randn(length(x),2));  %理论上dw_h + dw_v + dw_gama = w
        dw_poly_gama_fnc = mss2fnc(dw_poly_gama,x,randn(length(x),2));  %就是说三个偏导加起来就是全微分

        %% Put together，累加句柄
        W_exec = 'W_eval = @(ml)';  %函数句柄，W_eval将是m1操作的一个函数
        
        for i = 1:length(w_poly)
            if i<length(w_poly)
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d) +',i,i));      %水平串联字符串
            else
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d);',i,i));
            end
        end

        %% Execute，算好了保存输出
        eval(W_exec);   %computes the sum above，执行文本中的 MATLAB 表达式，文本转功能
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_h_fnc','dw_poly_v_fnc','dw_poly_gama_fnc','W_upper'); %都存进去，以供调用

    end
end
end
