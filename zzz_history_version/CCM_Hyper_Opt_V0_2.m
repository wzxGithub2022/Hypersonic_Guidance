%% �Ż��õ�CCM��function
% solved��һ����ʾSOS�����Ƿ�����flag
% �õĿ��ȫ��spotless��ֻ�ǵ�����Mosek���������
function [solved,w_lower,w_upper] = ...
    CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
    condn,lambda,ccm_eps,return_metric,varargin)
%% Current state load������״̬�����Ϣ
h1 = state_base(1);
v1 = state_base(2);
gama1 = state_base(3);
alpha1 = state_base(4);

global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

h0 = h1*R0;                      %����һ��ʾ��
v0 = v1*sqrt(R0*g);

%% State-space and dynamics��״̬�ռ�Ͷ���ѧ
W_scale = diag([0.001;0.001;0.001;0.001;0.001]);
norm_scale = 1e-4;

%handle function��������������������б�ѩ�����ʽ
sin_x = @(x) 0.985400513847416*x - 0.142496853019355*x^3;
cos_x = @(x) 0.999396553656190 - 0.495559134405119*x^2 + 0.036782872656059*x^4;

%states��״̬��y,h,v,gama,alpha
x = msspoly('x',5);         %ע�⣺��ЩxҲ�ǹ�һ����ģ��㶯��ѧϵ��ʱ��Ҫ����һ��

%pos_def indeterminates��������δȷ��ֵ
dfiv = msspoly('dfiv',5);   %d5��ϵͳ��5��״̬
dfor = msspoly('dfor',4);   %d4�����ƻ�Ӱ��һ��״̬

sin_g = sin_x(x(4));        
cos_g = cos_x(x(4));        %gama
sin_a = sin_x(x(5));
cos_a = cos_x(x(5));        %alpha

%dynamics f������ѧ
S   = 0.5026;                       %�ο����
% rou = 1.225 * exp(-x(2)/7110);    %�ܶ�rou��(msspoly���ܼ���ָ����̩��չ���ƽ�)
hf  = x(2)*R0/7110;                 %h fake���м����
rou = 1.225 * (1.124508748184077 + 0.459262515874491*hf + 3.155718757366057*hf^2 + 1.426742501850045*hf^3 + 0.374836249394692*hf^4);
q   = 0.5 * rou * (x(3)*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * x(3)*sqrt(R0*g);               %q fake��α��ѹ����Լ���ٶ�v
M   = x(3)*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                      %����

CL  = 0.4172 + 19.41*x(5) + 10.17*x(5)^2 - M*(0.1004 + 0.7536*x(5));
L_nor = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %����

%�б�ѩ��ƽ�1/V���Թ�һ��֮���V������[3,7]
V_division = 1.101028089323823 - 0.473922941895471*x(3) + 0.099720451270527*x(3)^2 - 0.010265506659430*x(3)^3 + 4.140771839924636e-04*x(3)^4;

f1 = -x(3) * cos_a;
f2 = x(3) * sin_a;
f3 = -D_nor - sin_g;
f4 = Lf_nor - cos_g*V_division;
f5 = 0;

f = [f1;f2;f3;f4;f5];               %����ѧ

%����΢�ֶ���ѧ��df_mat��ǰ���У����ƻ�Ӱ�쵽��һ�У������õ���
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
 
df_perp = [df1;df2;df3;df4];        %΢�ֶ���ѧ

B_perp = [eye(4);zeros(1,4)];       %�����Ǹ�B��ֱ

%% Initialize problem���������

prog = spotsosprog;
prog = prog.withIndeterminate(x);           %Indeterminate��������
prog = prog.withIndeterminate(dfiv);
prog = prog.withIndeterminate(dfor);

[prog, w_lower] = prog.newPos(1);           %positive��������һ������ֵ
[prog, w_upper] = prog.newPos(1);

%% Parametrize W����W������

w_states = x(2:4);                      %����˵w��Ϊһ��ccmӦ���Ǹ߶�h���ٶ�v���������gama��function������ʵ��CCM������֮һ������W����Ҫ��״̬��Ϊw_states
w_poly = monomials(w_states,0:6);       %����״̬����Ӧ��д��0�׵�6�׵ĵ���ʽ

W_list      = cell(length(w_poly),1);
W_perp_list = cell(length(w_poly),1);
W_pc_list   = cell(length(w_poly),1);
W_c_list    = cell(length(w_poly),1);

[prog, W_perp_list{1}] = prog.newSym(4);    %������4ά����ϵͳ���µ�4ά�Գƾ���symmetric
[prog, W_pc_list{1}] = prog.newFree(4,1);   %newSym���˻���newFree���µ����ɱ���
[prog, W_c_list{1}] = prog.newSym(1);
W_list{1} = [W_perp_list{1},W_pc_list{1};
             W_pc_list{1}', W_c_list{1}];   %��ô��ֻ

W = W_list{1}*w_poly(1);    %����һ��,����˵W��w_poly��ÿһ�������ϵ����w_poly��w_states��ÿһ���0��6�ε���ʽ������ϵ

for i = 2:length(w_poly)
    [prog, W_perp_list{i}] = prog.newSym(4);    %��ʽ�Ϻ͵�һ����ȫһ��
    [prog, W_pc_list{i}] = prog.newFree(4,1);
    [prog, W_c_list{i}] = prog.newSym(1);
    W_list{i} = [W_perp_list{i},W_pc_list{i};
        W_pc_list{i}', W_c_list{i}];
    
    W = W + W_list{i}*w_poly(i);    %���ۼ�
end

W_perp = W(1:4,1:4);                    %ȥ������ֱ��Ӱ���
dW_f_perp = diff(W_perp(:),x)*f;        %���diff����΢�֣�(:)�Ѿ�д��һ�е���ʽ��
dW_f_perp = reshape(dW_f_perp,4,4);     %�����˵��reshape�Ǹ��ö���

[prog, W_upper] = prog.newSym(n);       %W_upper��5*5msspoly���µĶԳƾ���

%% Definiteness conditions��Լ������
% �ⲿ�ֵ����Թ淶Ӧ�úܴ�̶��ϲο�mosek

%ע�⣺�����������CCM���Ƿ�Ӧ���������box����ÿ��״̬ƽ��㣬��h_lim^2 - (x(2)-x_standard)^2
%Lagrange multipliers���������ճ���
box_lim = [ h_lim^2 - (x(2)-h1)^2;
            v_lim^2 - (x(3)-v1)^2;
            gama_lim^2 - (x(4)-gama1)^2;
            alpha_lim^2 - (x(5)-alpha1)^2];       %�޶���״̬��Χ

l_order = 5;    %������޹ص�״̬���ף�dfiv��Ҫ2��
[pos_def_mon, pos_def_mat] = monomials([w_states;dfiv],0:l_order);      %monomials����ʽ��mat��¼�����ɹ���
pos_def_keep = find( sum( pos_def_mat(:, length(w_states)+1:end) ,2 ) == 2 );   %����dfiv��2��
pos_def_mon = pos_def_mon(pos_def_keep);            
[prog, Ll] = prog.newSOSPoly(pos_def_mon,3);
[prog, Lu] = prog.newSOSPoly(pos_def_mon,3);        %L1 Lu����2*1��msspoly��СL�ǵ�ֵ��u�Ǹ�ֵ���������µ�SOS��poly

lc_order = 6;   %����������ص�״̬���ף�dfiv��Ҫ2��
l_ccm_states = [x(2);x(3);x(4);x(5)];
[ccm_def_mon, ccm_def_mat] = monomials([l_ccm_states;dfor],0:lc_order);          
ccm_def_keep = find(sum(ccm_def_mat(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dfor
ccm_def_mon = ccm_def_mon(ccm_def_keep);
[prog, Lc] = prog.newSOSPoly(ccm_def_mon,4);

%W uniform bounds��W�ķ�ΧԼ��
prog = prog.withPos(w_lower-1);                     %w_lower > 1��������Լ��
prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper����������Լ��

%Condition bound��������Լ��
prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn��������Լ��

%W pos def��positive definite������Լ��
prog = prog.withSOS( (dfiv'*W*dfiv - w_lower*(dfiv'*dfiv)) - (Ll'*box_lim(1:3)) );      %Լ��ΪSOS
prog = prog.withSOS( dfiv'*(W_upper - W)*dfiv - (Lu'*box_lim(1:3)));                    %Լ��ΪSOS

%CCM condition��CCM��������֤����
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %������G��x���Ǹ�������ǰ�󶼳˸�B��ֱ��W�ͱ����W_perp
prog = prog.withSOS( (dfor'*R_CCM*dfor - ccm_eps*(dfor'*dfor)) - (Lc'*box_lim) );   %Լ��ΪSOS

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

        W_sol = zeros(n,n,length(w_poly));  
        NNZ_list = zeros(length(w_poly),1); 
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-3); %clean Remove terms with small coefficients
            if sum(sum(abs(W_sol(:,:,i)))) > 0      %����ֵ�����㣬֤����W_sol(:,:,i)��������ˣ�����
                NNZ_list(i) = 1;                    %non_negative��������һ�����
            end
        end
        w_poly = w_poly(find(NNZ_list));        %find���ص���������w_poly�ػ�����
        W_sol = W_sol(:,:,find(NNZ_list));      %W_solҲ�ػ�����������Ҳû�ų�ʲô
        
        fprintf('%d non-zero monomials\n',length(w_poly)); %ͳ�Ʒ��㣨��Сֵtolerance������
        
        dw_poly_h = diff(w_poly,x(2));   %�������Ӿ�����ƫ��
        dw_poly_v = diff(w_poly,x(3));
        dw_poly_gama = diff(w_poly,x(4));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-3);

        %% Create monomial functions������������
        w_poly_fnc       = mss2fnc(w_poly,    x,randn(length(x),2));  %computes all monomials
        dw_poly_h_fnc    = mss2fnc(dw_poly_h, x,randn(length(x),2));
        dw_poly_v_fnc    = mss2fnc(dw_poly_v, x,randn(length(x),2));  %������dw_h + dw_v + dw_gama = w
        dw_poly_gama_fnc = mss2fnc(dw_poly_gama,x,randn(length(x),2));  %����˵����ƫ������������ȫ΢��

        %% Put together���ۼӾ��
        W_exec = 'W_eval = @(ml)';  %���������W_eval����m1������һ������
        
        for i = 1:length(w_poly)
            if i<length(w_poly)
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d) +',i,i));      %ˮƽ�����ַ���
            else
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d);',i,i));
            end
        end

        %% Execute������˱������
        eval(W_exec);   %computes the sum above��ִ���ı��е� MATLAB ���ʽ���ı�ת����
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_h_fnc','dw_poly_v_fnc','dw_poly_gama_fnc','W_upper'); %�����ȥ���Թ�����

    end
end
end
