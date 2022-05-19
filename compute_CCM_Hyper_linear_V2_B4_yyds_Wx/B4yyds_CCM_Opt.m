%为了提高程序的可比性和运行效率，基于在线算法框架，对B4_fix进行重构
function [solved,w_lower,w_upper,W_upper] = B4yyds_CCM_Opt(n,A,b,condn,lambda,ccm_eps,return_metric,varargin)
%% State-space and dynamics，状态空间和动力学
W_scale = diag([0.1;0.1;0.01;0.01]);
norm_scale = 1e-4;  %这在优化里乘了sum(a)，权重是不是太大了

%states，状态：y,h,v,gama，4*1的msspoly
x = msspoly('x',4);         %注意：这些x也是归一化后的，算动力学系数时需要反归一化

%pos_def indeterminates，正定的未确定值
dfor = msspoly('dfor',4);   %d4，系统共4个状态
dtwo = msspoly('dtwo',2);   %d2，控制会影响2个状态

%----------------------BEGIN：动力学和微分动力学---------------------------%

%dynamics f，动力学
f = A * x;             %动力学，4*1，含x，这里的f得和x保持同样的维度

%基于微分动力学的df_mat的前2行，控制会影响到后2行，为了和B_perp配合
df_perp = A(1:2,:);        %微分动力学，2*4，含x，这里的df得和B_perp配合

%------------------------END：动力学和微分动力学---------------------------%

B_perp = [eye(2);zeros(2,2)];       %这是那个B垂直
%% Initialize problem，搭建问题框架

prog = spotsosprog;
prog = prog.withIndeterminate(x);           %Indeterminate，待定量
prog = prog.withIndeterminate(dfor);        %x是状态量
prog = prog.withIndeterminate(dtwo);        %d--是辅助变量

[prog, w_lower] = prog.newPos(1);           %positive，生成了一个正的值
[prog, w_upper] = prog.newPos(1);

%% Parametrize W，将W参数化
%之前w是高度h、速度v、弹道倾角gama的function，这样是不是削太多了
%不应该呀，你的射程状态应该是完全无关，或者说你写上也会自动变零，也没关系
w_states = x(1:2);                        %这其实是CCM的条件之一，生成W所需要的状态成为w_states
w_poly = monomials(w_states,0:4);       %只剩一个状态量，4阶甚至觉得高了

W_list      = cell(length(w_poly),1);
W_perp_list = cell(length(w_poly),1);
W_pc_list   = cell(length(w_poly),1);
W_c_list    = cell(length(w_poly),1);

[prog, W_perp_list{1}] = prog.newSym(2);    %定义了4维的新系统，新的4维对称矩阵symmetric
[prog, W_pc_list{1}] = prog.newFree(2,2);   %newSym是退化的newFree，新的自由变量
[prog, W_c_list{1}] = prog.newSym(2);
W_list{1} = [W_perp_list{1},W_pc_list{1};
             W_pc_list{1}', W_c_list{1}];   %这么大只

W = W_list{1}*w_poly(1);    %做第一个,就是说W和w_poly的每一项都建立关系，而w_poly和w_states的每一项的0到6次单项式建立关系

for i = 2:length(w_poly)
    [prog, W_perp_list{i}] = prog.newSym(2);    %形式上和第一个完全一样
    [prog, W_pc_list{i}] = prog.newFree(2,2);
    [prog, W_c_list{i}] = prog.newSym(2);
    W_list{i} = [W_perp_list{i},W_pc_list{i};
        W_pc_list{i}', W_c_list{i}];
    
    W = W + W_list{i}*w_poly(i);    %并累加
end

W_perp = W(1:2,1:2);                    %去除控制直接影响的
dW_f_perp = diff(W_perp(:),x)*f;        %这个diff是求微分，(:)已经写成一列的形式了，这里的f得和x保持同样的维度
dW_f_perp = reshape(dW_f_perp,2,2);     %吴恩达说过reshape是个好东西

[prog, W_upper] = prog.newSym(n);       %W_upper是4*4msspoly，新的对称矩阵

%% Definiteness conditions，约束条件
%-------------------------box_lim Start-----------------------------------%
% %不对，box_lim应该是一个大的范围，否则会有反向的风险，干脆删掉，写成全局的形式
% %注意：因为是线性化，所以误差小量化
% %Lagrange multipliers，拉格朗日乘子
% box_lim = [ h_lim^2 - x(2)^2;
%             v_lim^2 - x(3)^2;
%             gama_lim^2 - x(4)^2];       %限定了状态范围
% % %添加
% % prog = prog.withSOS(box_lim(1));
% % prog = prog.withSOS(box_lim(2));
% % prog = prog.withSOS(box_lim(3));
% %这里是约束W，而W只和w_states有关
% l_order = 4;    %与控制无关的状态2阶，dfor需要2阶
% [pos_def_mon, pos_def_mat] = monomials([w_states;dfor],0:l_order);      %monomials单项式，mat记录了生成规则
% pos_def_keep = find( sum( pos_def_mat(:, length(w_states)+1:end) ,2 ) == 2 );   %保留dfor的2阶
% pos_def_mon = pos_def_mon(pos_def_keep);            
% [prog, Ll] = prog.newSOSPoly(pos_def_mon,1);
% [prog, Lu] = prog.newSOSPoly(pos_def_mon,1);        %L1 Lu都是2*1的msspoly，小L是低值，u是高值，生成了新的SOS的poly
% %这里是约束G(x)/R_CCM，与影响系统的状态量都有关
% lc_order = 5;   %包含控制相关的状态3阶，dtwo需要2阶
% l_ccm_states = [x(2);x(3);x(4)];
% [ccm_def_mon, ccm_def_mat] = monomials([l_ccm_states;dtwo],0:lc_order);          
% ccm_def_keep = find(sum(ccm_def_mat(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dtwo
% ccm_def_mon = ccm_def_mon(ccm_def_keep);
% [prog, Lc] = prog.newSOSPoly(ccm_def_mon,3);
% %W uniform bounds，W的范围约束
% prog = prog.withPos(w_lower-1);                     %w_lower > 1，正条件约束，这个1也可以调整
% prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper，半正定条件约束
% %Condition bound，条件数约束
% prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn，正条件约束
% %W pos def，positive definite，正定约束，用SOS的方式来约束，在box范围内W都满足条件
% prog = prog.withSOS( (dfor'*W*dfor - w_lower*(dfor'*dfor)) - (Ll'*box_lim(1)) );      %约束 W > w_lower*eye4
% prog = prog.withSOS( dfor'*(W_upper - W)*dfor - (Lu'*box_lim(1)));                    %约束 W_upper > W
% %CCM condition，CCM条件，保证收缩
% R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %本质是G（x）那个条件，前后都乘个B垂直，W就变成了W_perp
% prog = prog.withSOS( (dtwo'*R_CCM*dtwo - ccm_eps*(dtwo'*dtwo)) - (Lc'*box_lim) );   %约束 R_CCM > epsilion*eye2
%-------------------------box_lim End-------------------------------------%

%-------------------------box_lim Start-----------------------------------%
%W uniform bounds，W的范围约束
prog = prog.withPos(w_lower-1);                     %w_lower > 1，正条件约束，这个1也可以调整
prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper，半正定条件约束
%Condition bound，条件数约束
prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn，正条件约束
%W pos def，positive definite，正定约束，用SOS的方式来约束，在box范围内W都满足条件
prog = prog.withSOS( dfor'*W*dfor - w_lower*(dfor'*dfor) );      %约束 W > w_lower*eye4
prog = prog.withSOS( dfor'*(W_upper - W)*dfor );                    %约束 W_upper > W
%CCM condition，CCM条件，保证收缩
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %本质是G（x）那个条件，前后都乘个B垂直，W就变成了W_perp
prog = prog.withSOS( dtwo'*R_CCM*dtwo - ccm_eps*(dtwo'*dtwo) );   %约束 R_CCM > epsilion*eye2
%-------------------------box_lim End-------------------------------------%

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
        
        clean_epsilon = 1e-6;

        W_sol = zeros(n,n,length(w_poly));  %从SOS_soln.eval(W_list{i})里取值并筛掉小量
        NNZ_list = zeros(length(w_poly),1); 
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),clean_epsilon); %clean Remove terms with small coefficients，1e-3是不是清理得太多了，OMG一动百动，control的问题在这里得到了解决
            if sum(sum(abs(W_sol(:,:,i)))) > 0      %绝对值大于零，证明该W_sol(:,:,i)算出东西了，保留
                NNZ_list(i) = 1;                    %non_negative，保留的一个标记
            end
        end
        w_poly = w_poly(find(NNZ_list));        %find返回的是索引，w_poly重获新生
        W_sol = W_sol(:,:,find(NNZ_list));      %W_sol也重获新生，好像也没排除什么
        
        fprintf('%d non-zero monomials\n',length(w_poly)); %统计非零（极小值tolerance）单项
        
        dw_poly_h = diff(w_poly,x(2));   %看这样子就是求偏导
        %寂寞标记
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),clean_epsilon);   %clean问题标注

        %% Create monomial functions，创造句柄函数
        w_poly_fnc       = mss2fnc(w_poly,    x,randn(length(x),2));  %computes all monomials
        dw_poly_h_fnc    = mss2fnc(dw_poly_h, x,randn(length(x),2));  %现在的w_poly只是h的函数了
        %寂寞标记

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
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_h_fnc','W_upper'); %都存进去，以供调用

    end
end

end