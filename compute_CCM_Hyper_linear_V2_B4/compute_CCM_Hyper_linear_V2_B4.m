%% compute_CCM_Hyper����
% �汾������˵����V2 B4
% box_limӦ��ɾ���ã���Ҫ��ɾ�����h(x)��û������
% ���Կ���ֹͣ״̬��չ�����������ù��ǣ����B��B_perp����Ҫ��д
% ddf*x��A*A���Թ�����Ϊ����

% clear all; close all; clc;
% yalmip('clear'); %��ͼ����Ĺ��ܣ����Ҳ���
% warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load������״̬�����Ϣ

load('CCM_Dynamitics_B.mat');
load('CCM_Dynamitics_df.mat');

global state_num state_CCM CCW_upper

% state_num = 1;  %������һ������Ҫ����һ�зų���

%���Ի�֮����Ϣ������df_full��B�У�state_base/f_ful/ddf_full��ûʲô��
i = state_num;
B_full = B_mat_value(:,:,i);        % 4*1
df_full = df_mat_value(:,:,i);      % 4*4 double

global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

%% Constants������

n = 4;                      %����ѧ��4����
ccm_eps = 0.05;             %���ֵ0.05�������ǲ���̫�Ͽ���

%ע�⣬����ÿһ��������
h_lim = 100/R0;                  
v_lim = 30/sqrt(R0*g);
gama_lim = pi/180;          %�Ƕȵĵ�λ�ǣ�����

%% Uncomment for global optimization of bound�������OPT-CCM

% Ѱ����С�������Ĺ��̣���1.2�����������ٶ��ַ�����
% �Ż����յ���ʲô���ڲ�ͬ���������£�CCM���������ʱ������ָ���ǲ�ͬ�ģ���һϵ��lambda��������Ӧ������ָ�꣬�˹���ѡ���ŵ�lambda
% eucǰ����sqrt(double(w_upper/w_lower))/lambda���ȼ��ڣ�1/lambda^2��*(w_upper/w_lower)��������d_bar * diag(W_upper)^1/2
% d_barsǰ����sqrt(double(1/w_lower))/lambda��������alpha_w/lambda

% lambda_range = linspace(0.1,1,10);              %���Բ�ֵ��lambda�ķ�Χ
lambda_range = 0.5;                               %��ס��
lambda_range = (1/100)*round(lambda_range*100);   %��������Բ��
euc_bounds = NaN(length(lambda_range),1);         %ע�⣬ǰ���Ż�Ѱ��CCM�ͺ���׼ȷ����ʱ�ĺ��岻�ϸ���ͬ
d_bars = NaN(length(lambda_range),1);             %bar������һ��
cond_bound = NaN(length(lambda_range),1);         %������

eps = 1;
condn_prev = 5;                                  %Ԥ��������
return_metric = 0;

for ll = 1:length(lambda_range)       %��ÿһ��lambda��������Сд�ġ�LL��,���ۣ�������
    
    lambda = lambda_range(ll);
    fprintf('**********\n');
    fprintf('lambda: %f\n', lambda);
    solved = 0;
    
    %Determine upper bound�������Ͻ�
    cond_l = condn_prev;              
    cond_u = 1.2*condn_prev;
%------------------------������һֱ��ѭ�����ǿ���������䣬sos����һֱ�ò������-----------------------%
    while (~solved) 
        fprintf(' cond_u: %.2f: ', cond_u);
        [sos_prob,~,~] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,...
                                df_full,cond_u,lambda,ccm_eps,return_metric);
        if (sos_prob == 0)
            solved = 1;
            fprintf('feasible \n');
        else
            %shift up condition number range
            %������˵Ļ���һ������condition number��ֻ����һ���������
            fprintf('�㿴��ⲻ�����˰�\n');
            cond_l = cond_u;
            cond_u = 1.2*cond_u;
        end
        if cond_l > 1000       %Ϊ�˲��ó������������������������ѭ����
            break;
        end
    end
%---------------------------------------------------------------------------------------------------%
    if (solved)           %ֱ������ˣ���euc_bounds(ll)������
        euc_bounds(ll) = sqrt(cond_u)/lambda;       %���������һ��ֵ����ʱ���޷�����d_bars
        fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
    else
        continue;
    end
    
    %Now do bisection search������������
    while(cond_u - cond_l >= eps)     %���ϴ�ϸ���Ż�
        condn = (cond_l+cond_u)/2;    %condnֻ�Ǹ��м����
        fprintf(' cond: %.2f', condn);
        [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,...
                                df_full,condn,lambda,ccm_eps,return_metric);
        if (sos_prob == 0)
            fprintf(' feasible\n');       %��Ϳ����ˣ�������ȡСһ�㣬���󶯵Ļ���euc��d�ͻ������Ż�һЩ      
            euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;   %euc�ȼ��ڣ�1/lambda^2��*(w_upper/w_lower)
            d_bars(ll) = sqrt(double(1/w_lower))/lambda;             %������״μ���d_bars
            cond_u = condn;
        else
            fprintf(' infeasible\n');     %�����еĻ������ҷŴ�һ��
            cond_l = condn;
        end
    end

%     condn_prev = cond_u;            %���д�������ķ��˼·��lambda����Ԥ��������Ҳ����
    cond_bound(ll) = cond_u;        %����С��ȫ�����

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end
% figure(state_num)
% subplot(211),plot(lambda_range,euc_bounds);title('Euc_bound');
% subplot(212),plot(lambda_range,d_bars);title('d_bars');
% pause; 

%% Pick a (lambda,condn)���������ȭ����һ�׵�

lambda = 0.5;          %������������
condn = cond_u;        %�̶������ʣ�����������Ӧ
return_metric = 1;

save_name = 'metric_Hyper_vectorized_linear.mat';          %֮��Ĳ���ֻ������һ��CCM_OTP������save_name�Ѷ���ȡ������֮ǰ����Ĳ����Ѿ�������Ż�CCM�Ĺ���
%%
[sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,...
                                df_full,condn,lambda,ccm_eps,return_metric,save_name);      %����һ����ã�sos_prob:sum of squares problem

%% Compute aux control bound���������Ʊ߽�
%�����뿴�����������(lambda��condn)����£�CCM��ʲô���ı��֣���Ҫ����u���ܿص�ס
%���ƺ��Ŷ�����ϵģ�����һ���Ŷ��Ͻ磬�㻹�����������ʲô��tube�ڣ����������u�����

load(save_name);
disp('Checking CCM conditions and Computing control bound...');

%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%

syms y h v gama           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ

%dynamics f������ѧ�����㿴͸�������f_mat�ͺ��ڵ�f(1...4)�Ĺ�ϵ����ȫ���ÿ�������compute��opt��д������
f_mat(y, h, v, gama) = df_full * [y; h; v; gama];               %����ѧ��4*1����x

%΢�ֶ���ѧ��df_mat
df44 = df_full * df_full;
df_mat = df44;            %4*4����΢��ϵ�����󣬲���delta_x

%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%

%%
%Bֻ��������Ͻ������ù������Կ����ӵ�alpha_dot�����Ǿ��ǿ���������ֱ�Ӹ�ǿ�����Ŀ��ƿ϶��Ǹ��ÿص�
B = B_full;       % 4*1 ��double

%дΪ0.1����������
load('CCM_Dynamitics_D.mat')
load('CCM_Dynamitics_Lf.mat')
D_full = D_mat_value(state_num);
Lf_full = Lf_mat_value(state_num);

Bw = [zeros(2,2);
      D_full*0.01, 0;
      0, Lf_full*0.01];
  
B_perp = [eye(2);zeros(2,2)];            %�����Ǹ�B��ֱ
%%
ctrl_N = 6;               %���Բ�ֵ�ĸ��������ǰ���ô�����״̬ȫģ����һ��
h_range = linspace(-h_lim, h_lim, ctrl_N);      %ע�⣬����״̬��ֵӦ��ÿ������״̬��������[h0-h_lim,h0+hlim]
v_range = linspace(-v_lim, v_lim, ctrl_N);
gama_range = linspace(-gama_lim, gama_lim, ctrl_N);

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N);       %��3ά���飬��ΪҪ����3��״̬
eig_CCM_min = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_CCM_max = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,2);                     %Wֻ��i������Ϊ��Ĺ��ǻ�Ӱ��v��gama����������Wֻ�Ǹ߶�h�ĺ���
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);
%%
for i = 1:length(h_range)
    fprintf(' i = %d , calculating \n',i);
    for j = 1:length(v_range)
        fprintf(' j = %d , calculating \n',j);
        for k = 1:length(gama_range)
                
                x = [randn(1,1);h_range(i);v_range(j);gama_range(k)];    %������������״̬����
                
                W = W_eval(w_poly_fnc(x));  %��״̬xִ��W�����ɹ���õ�W
                M = W\eye(n);               %����������
                Theta = chol(M);            %R = chol(A) ���ھ��� A �ĶԽ��ߺ������������������Ǿ��� R�����㷽�� R'*R=A��
                Theta_Bw = Theta*Bw;
                sigma_ThBw(i,j,k) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %Ϊ����ƫ���Ͻ��׼��
                
                L = chol(W);
                
                %��yû��ϵ��д�Ϸ���f��df��߳˳���Ҳ��0
                f = double( f_mat(0,x(2),x(3),x(4)) );           %״̬΢��
                df = df_mat;         %״̬����΢��
                
                % W�ǶԳƵģ�FҲ�ǶԳƵģ�������G��x����F��Ҫ����������Ķ�Ӧ��ϵ�й�bug���ȸ�Ϊ��
                F = -W_eval(dw_poly_h_fnc(x)) * f(2) + df*W + W*df' + 2*lambda*W;  %B4��ȫ΢�־�ֻ��һ��ƫ΢����
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));     %L��W�ķֽ�
                delta_u(i,j,k) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %���ǿ�����
                
                R_CCM = -B_perp'*F*B_perp;  %ǰ�󶼳�B��ֱ��R_CCM��Ҫ����
                
                eig_CCM_min(i,j,k) = min(eig(R_CCM));
                eig_CCM_max(i,j,k) = max(eig(R_CCM));
                eig_W(i,1) = min(eig(W));
                eig_W(i,2) = max(eig(W));
                
        end
    end
end
%%
alpha_w = max(sigma_ThBw(:));       %(:)��ʾ��������Ԫ�ز����б�ʾ
d_bar = alpha_w/lambda;             %��d_barû�г�w�Ͻ磬������J_CCMҲ�����ٳ�w�Ͻ��ˣ�һ����λ
disp('d_bar'); 
disp(d_bar);

disp('Control:'); 
disp(max(d_bar*delta_u(:)));        %���ǿ����Ͻ�

disp('W:'); 
disp(min( eig_W(:,1) ));      %С��С
disp(max( eig_W(:,2) ));      %���д�

disp('min eig CCM:'); 
disp(min(eig_CCM_min(:)));

disp('max eig CCM:'); 
disp(max(eig_CCM_max(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));
%% auto test all ����
state_CCM(1,state_num) = alpha_w/lambda;
state_CCM(2,state_num) = max(d_bar*delta_u(:));
state_CCM(3,state_num) = min(eig_W(:,1));
state_CCM(4,state_num) = max(eig_W(:,2));
state_CCM(5,state_num) = min(eig_CCM_min(:));
state_CCM(6,state_num) = max(eig_CCM_max(:));
state_CCM(7:10,state_num) = d_bar*sqrt(diag(W_upper));
CCW_upper(:,:,state_num) = W_upper;