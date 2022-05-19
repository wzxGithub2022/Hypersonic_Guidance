%% compute_CCM_Hyper����
% ����������ȫû���õ�spotless�Ķ���������ֻ���Ż��õĹ���
% ���⻹����opt��
% ����˵�����Ի�֮���൱��ϵͳ�����ˣ��˰汾����
clear all; close all; clc;
% yalmip('clear'); %��ͼ����Ĺ��ܣ����Ҳ���
warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load������״̬�����Ϣ

load('Trajectory_normalization.mat');
load('CCM_Dynamitics_f.mat');
load('CCM_Dynamitics_df.mat');
load('CCM_Dynamitics_ddf.mat');

%���Ի�֮����Ϣ������df_full��ddf_full�У�state_base��f_full��ûʲô��
i = 1;
h1 = Trajectory_normalization(i,3);
v1 = Trajectory_normalization(i,4);
gama1 = Trajectory_normalization(i,5);
alpha1 = Trajectory_normalization(i,6);
state_base = [h1,v1,gama1,alpha1];

f_full = f_mat_value(:,:,i);        % 5*1 double�����ֻ�Ƿ������ûʲô����
df_full = df_mat_value(:,:,i);      % 5*5 double
ddf_full = ddf_mat_value(:,:,:,i);     % 5*5*5

global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

%% Constants������

n = 5;                      %����ѧ��5����
ccm_eps = 0.05;             %���ֵ0.05�������ǲ���̫�Ͽ���

%ע�⣬����ÿһ��������
h_lim = 100/R0;                  
v_lim = 30/sqrt(R0*g);
gama_lim = pi/180;           %�������Ƕȵĵ�λ���ǣ�����
alpha_lim = pi/360;

%% Uncomment for global optimization of bound�������OPT-CCM

%Ѱ����С�������Ĺ��̣���1.2�����������ٶ��ַ�����
%�Ż����յ���ʲô���ڲ�ͬ���������£�CCM���������ʱ������ָ���ǲ�ͬ�ģ���һϵ��lambda��������Ӧ������ָ�꣬�˹���ѡ���ŵ�lambda
%eucǰ����sqrt(double(w_upper/w_lower))/lambda���ȼ��ڣ�1/lambda^2��*(w_upper/w_lower)��������d_bar * diag(W_upper)^1/2
%d_barsǰ����sqrt(double(1/w_lower))/lambda��������alpha_w/lambda

lambda_range = linspace(0.5,1,10);              %���Բ�ֵ��lambda�ķ�Χ
lambda_range = (1/100)*round(lambda_range*100);   %��������Ϊ����
euc_bounds = NaN(length(lambda_range),1);         %ע�⣬ǰ���Ż�Ѱ��CCM�ͺ���׼ȷ����ʱ�ĺ��岻�ϸ���ͬ��
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
        [sos_prob,~,~] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                f_full,df_full,ddf_full,cond_u,lambda,ccm_eps,return_metric);
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
        if cond_l > 100
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
        condn = (cond_l+cond_u)/2;
        fprintf(' cond: %.2f', condn);
        [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                f_full,df_full,ddf_full,condn,lambda,ccm_eps,return_metric);
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

    condn_prev = cond_u;              %����С��ȫ�����
    cond_bound(ll) = cond_u;

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end

pause; 

%% Pick a (lambda,condn)���������ȭ����һ�׵�

lambda = 0.7;          %������������
condn = 10;
return_metric = 1;

save_name = 'metric_Hyper_vectorized_linear.mat';          %֮��Ĳ���ֻ������һ��CCM_OTP������save_name�Ѷ���ȡ������֮ǰ����Ĳ����Ѿ�������Ż�CCM�Ĺ���
%%
[sos_prob, w_lower, w_upper] = CCM_Hyper_Opt_linear(n,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                f_full,df_full,ddf_full,condn,lambda,ccm_eps,return_metric,save_name);      %����һ����ã�sos_prob:sum of squares problem

%% Compute aux control bound���������Ʊ߽�
%�����뿴�����������(lambda��condn)����£�CCM��ʲô���ı��֣���Ҫ����u���ܿص�ס
%���ƺ��Ŷ�����ϵģ�����һ���Ŷ��Ͻ磬�㻹�����������ʲô��tube�ڣ����������u�����

load(save_name);
disp('Checking CCM conditions and Computing control bound...');

B = [zeros(4,1);1];       % 5*1 ��double
 
Bw = @(x)[zeros(2,2);     % �������Ŷ�����˼�ǣ��������������Ϊ0.01��������
           0.01/1600,0;
     0,0.01/(1600*2000);
                0,0];      %������д��ֱ�Ӿ���Сfunction�����������Ŷ���ϵͳ�����룬��СΪ5*2
  
B_perp = [eye(4);zeros(1,4)];            %�����Ǹ�B��ֱ

ctrl_N = 6;               %���Բ�ֵ�ĸ��������ǰ���ô�����״̬ȫģ����һ��
h_range = linspace(-h_lim, h_lim, ctrl_N);      %ע�⣬����״̬��ֵӦ��ÿ������״̬��������[h0-h_lim,h0+hlim]
v_range = linspace(-v_lim, v_lim, ctrl_N);
gama_range = linspace(-gama_lim, gama_lim, ctrl_N);
alpha_range = linspace(-alpha_lim, alpha_lim, ctrl_N);

%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%

syms h v gama alpha           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ

% %dynamics f������ѧ

f_mat(h, v, gama, alpha) = df_full(2:5,2:5) * [h; v; gama; alpha];               %����ѧ��4*1����x

%����΢�ֶ���ѧ��df_mat��ǰ���У����ƻ�Ӱ�쵽��һ�У������õ���
df1 = ddf_full(:,:,1) * [0; h; v; gama; alpha];
df2 = ddf_full(:,:,2) * [0; h; v; gama; alpha];
df3 = ddf_full(:,:,3) * [0; h; v; gama; alpha];
df4 = ddf_full(:,:,4) * [0; h; v; gama; alpha];

df1 = reshape(df1,[1,5]);
df2 = reshape(df2,[1,5]);
df3 = reshape(df3,[1,5]);
df4 = reshape(df4,[1,5]);
df5 = [0,0,0,0,0];

df_mat(h, v, gama, alpha) = [df1;df2;df3;df4;df5];            %5*5����x

%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);       %����ά���飬��Ϊ��4��״̬
eig_CCM_min = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_CCM_max = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,2);                     %Wֻ��i,j,k������Ϊ���u��Ӱ�칥�ǣ���������Wֻ�Ǹ߶ȡ��ٶȡ�������ǵĺ���
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);
%%
for i = 1:length(h_range)
    fprintf(' i = %d , calculating \n',i);
    for j = 1:length(v_range)
        fprintf(' j = %d , calculating \n',j);
        for k = 1:length(gama_range)
            for l = 1:length(alpha_range)      %����Сд��L��
                
                x = [randn(1,1);h_range(i);v_range(j);gama_range(k);alpha_range(l)];    %������������״̬����
                
                W = W_eval(w_poly_fnc(x));
                M = W\eye(n);               %����������
                Theta = chol(M);            %R = chol(A) ���ھ��� A �ĶԽ��ߺ������������������Ǿ��� R�����㷽�� R'*R=A��
                Theta_Bw = Theta*Bw(x);
                sigma_ThBw(i,j,k,l) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %Ϊ����ƫ���Ͻ��׼��
                
                L_nor = chol(W);
                
                f = double( f_mat(x(2),x(3),x(4),x(5)) );           %״̬΢��
                df = double( df_mat(x(2),x(3),x(4),x(5)) );         %״̬����΢��
                
                % W�ǶԳƵģ�FҲ�ǶԳƵģ�������G��x����F��Ҫ�������ȸ�Ϊ��
                F = -W_eval(dw_poly_h_fnc(x)) * f(1) - W_eval(dw_poly_v_fnc(x)) * f(2) - W_eval(dw_poly_gama_fnc(x)) * f(3)...      
                            + df*W + W*df' + 2*lambda*W;  %����ƫ΢���Ѿ�������ȫ΢�֣�alpha��ȫ��0��ֱ�Ӳ���д
                
                delta_u_den = eig((inv(L_nor))'*(B*B')*inv(L_nor));     %L��W�ķֽ�
                delta_u(i,j,k,l) = 0.5*max(eig((inv(L_nor))'*F*inv(L_nor)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %���ǿ�����
                
                R_CCM = -B_perp'*F*B_perp;  %ǰ�󶼳�B��ֱ��R_CCM��Ҫ����
                
                eig_CCM_min(i,j,k,l) = min(eig(R_CCM));
                eig_CCM_max(i,j,k,l) = max(eig(R_CCM));
                eig_W(i,j,k,1) = min(eig(W));
                eig_W(i,j,k,2) = max(eig(W));
                
            end
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
disp(min(min(min(eig_W(:,:,:,1)))));      %С��С
disp(min(max(max(eig_W(:,:,:,2)))));      %���д�

disp('min eig CCM:'); 
disp(min(eig_CCM_min(:)));

disp('max eig CCM:'); 
disp(max(eig_CCM_max(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));
