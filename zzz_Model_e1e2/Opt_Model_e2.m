%% ��ģ��e2
function [solved,w_lower,w_upper] = Opt_Model_e2(n,state_lim,condn,lambda,ccm_eps,return_metric,varargin)

%% Current state load������״̬�����Ϣ
global R0 g

%% State-space and dynamics��״̬�ռ�Ͷ���ѧ

W_scale = diag([0.001;0.001;0.001;0.001]);
norm_scale = 1e-4;

%states��״̬��y,h,v,gama,alpha
x = msspoly('x',4);         %ע�⣺��ЩxҲ�ǹ�һ����ģ��㶯��ѧϵ��ʱ��Ҫ����һ��

%pos_def indeterminates��������δȷ��ֵ
dfor = msspoly('dfor',4);   %d4��ϵͳ��4��״̬
dthr = msspoly('dthr',3);   %d3�����ƻ�Ӱ��һ��״̬

%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%

%dynamics f������ѧ
m  = 1600;        %����

CL = 0.4 * 180 / pi;
Cd = 0.02988 * CL^2;

%handle function��������������������б�ѩ�����ʽ
sin_x = @(x) 0.985400513847416*x - 0.142496853019355*x^3;
cos_x = @(x) 0.999396553656190 - 0.495559134405119*x^2 + 0.036782872656059*x^4;

sin_g = sin_x(x(4));        
cos_g = cos_x(x(4));        %gama

%�б�ѩ��ƽ�1/V���Թ�һ��֮���V������[3,7]
V_division = 1.101028089323823 - 0.473922941895471*x(3) + 0.099720451270527*x(3)^2 - 0.010265506659430*x(3)^3 + 4.140771839924636e-04*x(3)^4;

f1 = -x(3) * cos_g;
f2 = x(3) * sin_g;
f3 = - Cd/(m*g) - sin_g;
f4 = - cos_g * V_division;  %����ȫ�ֵ�˼·��alphaд������

f = [f1;f2;f3;f4];               %����ѧ�������f�ú�x����ͬ����ά��

%����΢�ֶ���ѧ��df_mat��ǰ���У����ƻ�Ӱ�쵽��һ�У������õ���
df11 = 0;
df12 = 0;
df13 = -cos_g;
df14 = x(3) * sin_g;
 df1 = [df11,df12,df13,df14];

df21 = 0;
df22 = 0;
df23 = sin_g;
df24 = x(3) * cos_g;
 df2 = [df21,df22,df23,df24];

df31 = 0;
df32 = 0;
df33 = 0;
df34 = -cos_g;
 df3 = [df31,df32,df33,df34];
 
df_perp = [df1;df2;df3];        %΢�ֶ���ѧ�������df�ú�B_perp���

%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%

B_perp = [eye(3);zeros(1,3)];       %�����Ǹ�B��ֱ

%% Initialize problem���������

prog = spotsosprog;
prog = prog.withIndeterminate(x);           %Indeterminate��������
prog = prog.withIndeterminate(dfor);        %x��״̬��
prog = prog.withIndeterminate(dthr);        %d4 d3�Ǹ�������

[prog, w_lower] = prog.newPos(1);           %positive��������һ������ֵ
[prog, w_upper] = prog.newPos(1);

%% Parametrize W����W������

w_states = x(3);                      %����˵w��Ϊһ��ccmӦ�����ٶ�v���������gama��function������ʵ��CCM������֮һ������W����Ҫ��״̬��Ϊw_states
w_poly = monomials(w_states,0:4);       %����״̬����Ӧ��д��0�׵�6�׵ĵ���ʽ��0:6--84��0:4--35

W_list      = cell(length(w_poly),1);
W_perp_list = cell(length(w_poly),1);
W_pc_list   = cell(length(w_poly),1);
W_c_list    = cell(length(w_poly),1);

[prog, W_perp_list{1}] = prog.newSym(3);    %������4ά����ϵͳ���µ�4ά�Գƾ���symmetric
[prog, W_pc_list{1}] = prog.newFree(3,1);   %newSym���˻���newFree���µ����ɱ���
[prog, W_c_list{1}] = prog.newSym(1);
W_list{1} = [W_perp_list{1},W_pc_list{1};
             W_pc_list{1}', W_c_list{1}];   %��ô��ֻ

W = W_list{1}*w_poly(1);    %����һ��,����˵W��w_poly��ÿһ�������ϵ����w_poly��w_states��ÿһ���0��6�ε���ʽ������ϵ

for i = 2:length(w_poly)
    [prog, W_perp_list{i}] = prog.newSym(3);    %��ʽ�Ϻ͵�һ����ȫһ��
    [prog, W_pc_list{i}] = prog.newFree(3,1);
    [prog, W_c_list{i}] = prog.newSym(1);
    W_list{i} = [W_perp_list{i},W_pc_list{i};
        W_pc_list{i}', W_c_list{i}];
    
    W = W + W_list{i}*w_poly(i);    %���ۼ�
end

W_perp = W(1:3,1:3);                    %ȥ������ֱ��Ӱ���
dW_f_perp = diff(W_perp(:),x)*f;        %���diff����΢�֣�(:)�Ѿ�д��һ�е���ʽ�ˣ������f�ú�x����ͬ����ά��
dW_f_perp = reshape(dW_f_perp,3,3);     %�����˵��reshape�Ǹ��ö���

[prog, W_upper] = prog.newSym(n);       %W_upper��4*4msspoly���µĶԳƾ���

%% Definiteness conditions��Լ������
%---------------------------------keep box--------------------------------%
v_lim_l = state_lim(1);
v_lim_u = state_lim(2);
v_sum = (v_lim_l + v_lim_u) / 2;
v_minus = (v_lim_u - v_lim_l) / 2;

gama_lim_l = state_lim(3);           %�������Ƕȵĵ�λ���ǣ�����
gama_lim_u = state_lim(4);
gama_sum = (gama_lim_l + gama_lim_u) / 2;
gama_minus = (gama_lim_u - gama_lim_l) / 2;

%Lagrange multipliers���������ճ���
box_lim = [ v_minus^2 - (x(3)-v_sum)^2;
            gama_minus^2 - (x(4)-gama_sum)^2];       %�޶���״̬��Χ
        
% prog = prog.withSOS( box_lim(1) );   
% prog = prog.withSOS( box_lim(2) ); 

%������Լ��W����Wֻ��w_states�й�
l_order = 3;    %������޹ص�״̬1�ף�dfor��Ҫ2��
[pos_def_mon, pos_def_mat] = monomials([w_states;dfor],0:l_order);      %monomials����ʽ��mat��¼�����ɹ���
pos_def_keep = find( sum( pos_def_mat(:, length(w_states)+1:end) ,2 ) == 2 );   %����dfiv��2��
pos_def_mon = pos_def_mon(pos_def_keep);            
[prog, Ll] = prog.newSOSPoly(pos_def_mon,1);
[prog, Lu] = prog.newSOSPoly(pos_def_mon,1);        %L1 Lu����2*1��msspoly��СL�ǵ�ֵ��u�Ǹ�ֵ���������µ�SOS��poly

%������Լ��G(x)/R_CCM����Ӱ��ϵͳ��״̬�����й�
lc_order = 4;   %����������ص�״̬2�ף�dfor��Ҫ2��
l_ccm_states = [x(3);x(4)];
[ccm_def_mon, ccm_def_mat] = monomials([l_ccm_states;dthr],0:lc_order);          
ccm_def_keep = find(sum(ccm_def_mat(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dfor
ccm_def_mon = ccm_def_mon(ccm_def_keep);
[prog, Lc] = prog.newSOSPoly(ccm_def_mon,2);

%W uniform bounds��W�ķ�ΧԼ��
prog = prog.withPos(w_lower-1);                     %w_lower > 1��������Լ��
prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper����������Լ��

%Condition bound��������Լ��
prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn��������Լ��

%W pos def��positive definite������Լ������SOS�ķ�ʽ��Լ������box��Χ��W����������
prog = prog.withSOS( (dfor'*W*dfor - w_lower*(dfor'*dfor)) - (Ll'*box_lim(1)) );      %Լ�� W > w_lower*eye5
prog = prog.withSOS( dfor'*(W_upper - W)*dfor - (Lu'*box_lim(1)));                    %Լ�� W_upper > W

%CCM condition��CCM��������֤����
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %������G��x���Ǹ�������ǰ�󶼳˸�B��ֱ��W�ͱ����W_perp
prog = prog.withSOS( (dthr'*R_CCM*dthr - ccm_eps*(dthr'*dthr)) - (Lc'*box_lim) );   %Լ�� R_CCM > epsilion*eye4
%---------------------------------keep box--------------------------------%

%--------------------------------throw box--------------------------------%
% %W uniform bounds��W�ķ�ΧԼ��
% prog = prog.withPos(w_lower-1);                     %w_lower > 1��������Լ��
% prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper����������Լ��
% %Condition bound��������Լ��
% prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn��������Լ��
% %W pos def��positive definite������Լ������SOS�ķ�ʽ��Լ������box��Χ��W����������
% prog = prog.withSOS( (dfor'*W*dfor - w_lower*(dfor'*dfor)) );      %Լ�� W > w_lower*eye5
% prog = prog.withSOS( dfor'*(W_upper - W)*dfor );                    %Լ�� W_upper > W
% %CCM condition��CCM��������֤����
% R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %������G��x���Ǹ�������ǰ�󶼳˸�B��ֱ��W�ͱ����W_perp
% prog = prog.withSOS( (dthr'*R_CCM*dthr - ccm_eps*(dthr'*dthr)) );   %Լ�� R_CCM > epsilion*eye4
%--------------------------------throw box--------------------------------%

%Լ����ɣ�����SDP�붨�滮����
options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint������Լ����ѹ����ķ�Χ��ʹ����������
free_vars = [prog.coneVar(2:end); prog.freeVar];   
len = length(free_vars);
[prog, a] = prog.newPos(len);       %a��free_vars������ͬ��ά�ȣ�������Լ��
prog = prog.withPos(-free_vars + a);    
prog = prog.withPos(free_vars + a);

%try-catch-end��ִ����䲢��������Ĵ���
    SOS_soln = prog.minimize(trace(W_scale*W_upper) + norm_scale*sum(a), @spot_mosek, options);     %����mosek��ֱ��һ��minimize�Ż�������
    %����һ���Ļ���������Ż�Ӧ���ܲ�������������ϰ���������Ϊ������Ʋ���
try
    %SOS_soln����prog
    solved = ~SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');    %�Ƚ��ַ��������Ƚϳɹ�ʱ��solved����0
    w_lower = double(SOS_soln.eval(w_lower));
    w_upper = double(SOS_soln.eval(w_upper));
catch
    %failed������ζ��CCMʧ�ܣ����ܽ������
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    return;
end

%% Parse���������ѽ����ȡ����

W_upper_mat = zeros(n);     %��ʵû���õ��ϣ���3D-Quadrotor�õõ�

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');

        W_sol = zeros(n,n,length(w_poly));  %��SOS_soln.eval(W_list{i})��ȡֵ��ɸ��С��
        NNZ_list = zeros(length(w_poly),1); 
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-9); %clean Remove terms with small coefficients
            if sum(sum(abs(W_sol(:,:,i)))) > 0      %����ֵ�����㣬֤����W_sol(:,:,i)��������ˣ�����
                NNZ_list(i) = 1;                    %non_negative��������һ�����
            end
        end
        w_poly = w_poly(find(NNZ_list));        %find���ص���������w_poly�ػ�����
        W_sol = W_sol(:,:,find(NNZ_list));      %W_solҲ�ػ�����������Ҳû�ų�ʲô
        
        fprintf('%d non-zero monomials\n',length(w_poly)); %ͳ�Ʒ��㣨��Сֵtolerance������
        
        dw_poly_v = diff(w_poly,x(3));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-9);

        %% Create monomial functions������������
        w_poly_fnc       = mss2fnc(w_poly,    x,randn(length(x),2));  %computes all monomials
        dw_poly_v_fnc    = mss2fnc(dw_poly_v, x,randn(length(x),2));  %������dw_h + dw_v + dw_gama = w

        %% Put together���ۼӾ��
        W_exec = 'W_eval = @(ml)';  %W_eval��ִ�ж���(m1)������������ƫ��
        
        for i = 1:length(w_poly)
            if i<length(w_poly)
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d) +',i,i));      %ˮƽ�����ַ���
            else
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d);',i,i));
            end
        end

        %% Execute������˱������
        eval(W_exec);   %computes the sum above��ִ���ı��е� MATLAB ���ʽ���ı�ת����
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_v_fnc','W_upper'); %�����ȥ���Թ�����

    end
end
end
