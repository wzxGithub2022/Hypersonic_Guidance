%Ϊ����߳���Ŀɱ��Ժ�����Ч�ʣ����������㷨��ܣ���B4_fix�����ع�
function [solved,w_lower,w_upper,W_upper] = B4yyds_CCM_Opt(n,A,b,condn,lambda,ccm_eps,return_metric,varargin)
%% State-space and dynamics��״̬�ռ�Ͷ���ѧ
W_scale = diag([0.1;0.1;0.01;0.01]);
norm_scale = 1e-4;  %�����Ż������sum(a)��Ȩ���ǲ���̫����

%states��״̬��y,h,v,gama��4*1��msspoly
x = msspoly('x',4);         %ע�⣺��ЩxҲ�ǹ�һ����ģ��㶯��ѧϵ��ʱ��Ҫ����һ��

%pos_def indeterminates��������δȷ��ֵ
dfor = msspoly('dfor',4);   %d4��ϵͳ��4��״̬
dtwo = msspoly('dtwo',2);   %d2�����ƻ�Ӱ��2��״̬

%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%

%dynamics f������ѧ
f = A * x;             %����ѧ��4*1����x�������f�ú�x����ͬ����ά��

%����΢�ֶ���ѧ��df_mat��ǰ2�У����ƻ�Ӱ�쵽��2�У�Ϊ�˺�B_perp���
df_perp = A(1:2,:);        %΢�ֶ���ѧ��2*4����x�������df�ú�B_perp���

%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%

B_perp = [eye(2);zeros(2,2)];       %�����Ǹ�B��ֱ
%% Initialize problem���������

prog = spotsosprog;
prog = prog.withIndeterminate(x);           %Indeterminate��������
prog = prog.withIndeterminate(dfor);        %x��״̬��
prog = prog.withIndeterminate(dtwo);        %d--�Ǹ�������

[prog, w_lower] = prog.newPos(1);           %positive��������һ������ֵ
[prog, w_upper] = prog.newPos(1);

%% Parametrize W����W������
%֮ǰw�Ǹ߶�h���ٶ�v���������gama��function�������ǲ�����̫����
%��Ӧ��ѽ��������״̬Ӧ������ȫ�޹أ�����˵��д��Ҳ���Զ����㣬Ҳû��ϵ
w_states = x(1:2);                        %����ʵ��CCM������֮һ������W����Ҫ��״̬��Ϊw_states
w_poly = monomials(w_states,0:4);       %ֻʣһ��״̬����4���������ø���

W_list      = cell(length(w_poly),1);
W_perp_list = cell(length(w_poly),1);
W_pc_list   = cell(length(w_poly),1);
W_c_list    = cell(length(w_poly),1);

[prog, W_perp_list{1}] = prog.newSym(2);    %������4ά����ϵͳ���µ�4ά�Գƾ���symmetric
[prog, W_pc_list{1}] = prog.newFree(2,2);   %newSym���˻���newFree���µ����ɱ���
[prog, W_c_list{1}] = prog.newSym(2);
W_list{1} = [W_perp_list{1},W_pc_list{1};
             W_pc_list{1}', W_c_list{1}];   %��ô��ֻ

W = W_list{1}*w_poly(1);    %����һ��,����˵W��w_poly��ÿһ�������ϵ����w_poly��w_states��ÿһ���0��6�ε���ʽ������ϵ

for i = 2:length(w_poly)
    [prog, W_perp_list{i}] = prog.newSym(2);    %��ʽ�Ϻ͵�һ����ȫһ��
    [prog, W_pc_list{i}] = prog.newFree(2,2);
    [prog, W_c_list{i}] = prog.newSym(2);
    W_list{i} = [W_perp_list{i},W_pc_list{i};
        W_pc_list{i}', W_c_list{i}];
    
    W = W + W_list{i}*w_poly(i);    %���ۼ�
end

W_perp = W(1:2,1:2);                    %ȥ������ֱ��Ӱ���
dW_f_perp = diff(W_perp(:),x)*f;        %���diff����΢�֣�(:)�Ѿ�д��һ�е���ʽ�ˣ������f�ú�x����ͬ����ά��
dW_f_perp = reshape(dW_f_perp,2,2);     %�����˵��reshape�Ǹ��ö���

[prog, W_upper] = prog.newSym(n);       %W_upper��4*4msspoly���µĶԳƾ���

%% Definiteness conditions��Լ������
%-------------------------box_lim Start-----------------------------------%
% %���ԣ�box_limӦ����һ����ķ�Χ��������з���ķ��գ��ɴ�ɾ����д��ȫ�ֵ���ʽ
% %ע�⣺��Ϊ�����Ի����������С����
% %Lagrange multipliers���������ճ���
% box_lim = [ h_lim^2 - x(2)^2;
%             v_lim^2 - x(3)^2;
%             gama_lim^2 - x(4)^2];       %�޶���״̬��Χ
% % %���
% % prog = prog.withSOS(box_lim(1));
% % prog = prog.withSOS(box_lim(2));
% % prog = prog.withSOS(box_lim(3));
% %������Լ��W����Wֻ��w_states�й�
% l_order = 4;    %������޹ص�״̬2�ף�dfor��Ҫ2��
% [pos_def_mon, pos_def_mat] = monomials([w_states;dfor],0:l_order);      %monomials����ʽ��mat��¼�����ɹ���
% pos_def_keep = find( sum( pos_def_mat(:, length(w_states)+1:end) ,2 ) == 2 );   %����dfor��2��
% pos_def_mon = pos_def_mon(pos_def_keep);            
% [prog, Ll] = prog.newSOSPoly(pos_def_mon,1);
% [prog, Lu] = prog.newSOSPoly(pos_def_mon,1);        %L1 Lu����2*1��msspoly��СL�ǵ�ֵ��u�Ǹ�ֵ���������µ�SOS��poly
% %������Լ��G(x)/R_CCM����Ӱ��ϵͳ��״̬�����й�
% lc_order = 5;   %����������ص�״̬3�ף�dtwo��Ҫ2��
% l_ccm_states = [x(2);x(3);x(4)];
% [ccm_def_mon, ccm_def_mat] = monomials([l_ccm_states;dtwo],0:lc_order);          
% ccm_def_keep = find(sum(ccm_def_mat(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dtwo
% ccm_def_mon = ccm_def_mon(ccm_def_keep);
% [prog, Lc] = prog.newSOSPoly(ccm_def_mon,3);
% %W uniform bounds��W�ķ�ΧԼ��
% prog = prog.withPos(w_lower-1);                     %w_lower > 1��������Լ�������1Ҳ���Ե���
% prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper������������Լ��
% %Condition bound��������Լ��
% prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn��������Լ��
% %W pos def��positive definite������Լ������SOS�ķ�ʽ��Լ������box��Χ��W����������
% prog = prog.withSOS( (dfor'*W*dfor - w_lower*(dfor'*dfor)) - (Ll'*box_lim(1)) );      %Լ�� W > w_lower*eye4
% prog = prog.withSOS( dfor'*(W_upper - W)*dfor - (Lu'*box_lim(1)));                    %Լ�� W_upper > W
% %CCM condition��CCM��������֤����
% R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %������G��x���Ǹ�������ǰ�󶼳˸�B��ֱ��W�ͱ����W_perp
% prog = prog.withSOS( (dtwo'*R_CCM*dtwo - ccm_eps*(dtwo'*dtwo)) - (Lc'*box_lim) );   %Լ�� R_CCM > epsilion*eye2
%-------------------------box_lim End-------------------------------------%

%-------------------------box_lim Start-----------------------------------%
%W uniform bounds��W�ķ�ΧԼ��
prog = prog.withPos(w_lower-1);                     %w_lower > 1��������Լ�������1Ҳ���Ե���
prog = prog.withPSD(w_upper*eye(n)-W_upper);        %w_upper*eye(n) > W_upper������������Լ��
%Condition bound��������Լ��
prog = prog.withPos(condn*w_lower - w_upper);       %w_upper/w_lower>condn��������Լ��
%W pos def��positive definite������Լ������SOS�ķ�ʽ��Լ������box��Χ��W����������
prog = prog.withSOS( dfor'*W*dfor - w_lower*(dfor'*dfor) );      %Լ�� W > w_lower*eye4
prog = prog.withSOS( dfor'*(W_upper - W)*dfor );                    %Լ�� W_upper > W
%CCM condition��CCM��������֤����
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);    %������G��x���Ǹ�������ǰ�󶼳˸�B��ֱ��W�ͱ����W_perp
prog = prog.withSOS( dtwo'*R_CCM*dtwo - ccm_eps*(dtwo'*dtwo) );   %Լ�� R_CCM > epsilion*eye2
%-------------------------box_lim End-------------------------------------%

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
        
        clean_epsilon = 1e-6;

        W_sol = zeros(n,n,length(w_poly));  %��SOS_soln.eval(W_list{i})��ȡֵ��ɸ��С��
        NNZ_list = zeros(length(w_poly),1); 
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),clean_epsilon); %clean Remove terms with small coefficients��1e-3�ǲ��������̫���ˣ�OMGһ���ٶ���control������������õ��˽��
            if sum(sum(abs(W_sol(:,:,i)))) > 0      %����ֵ�����㣬֤����W_sol(:,:,i)��������ˣ�����
                NNZ_list(i) = 1;                    %non_negative��������һ�����
            end
        end
        w_poly = w_poly(find(NNZ_list));        %find���ص���������w_poly�ػ�����
        W_sol = W_sol(:,:,find(NNZ_list));      %W_solҲ�ػ�����������Ҳû�ų�ʲô
        
        fprintf('%d non-zero monomials\n',length(w_poly)); %ͳ�Ʒ��㣨��Сֵtolerance������
        
        dw_poly_h = diff(w_poly,x(2));   %�������Ӿ�����ƫ��
        %��į���
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),clean_epsilon);   %clean�����ע

        %% Create monomial functions������������
        w_poly_fnc       = mss2fnc(w_poly,    x,randn(length(x),2));  %computes all monomials
        dw_poly_h_fnc    = mss2fnc(dw_poly_h, x,randn(length(x),2));  %���ڵ�w_polyֻ��h�ĺ�����
        %��į���

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
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_h_fnc','W_upper'); %�����ȥ���Թ�����

    end
end

end